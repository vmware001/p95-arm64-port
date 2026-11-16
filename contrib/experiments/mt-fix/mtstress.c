/*----------------------------------------------------------------------
| mtstress.c - differential multithreaded-gwnum stress test
|
| Purpose: find the intermittent WRONG-RESULT bug in the gwnum ARM64
| "asm converted to C" port when a single worker uses several threads
| (mprime CoresPerWorker > 1).  The bug shows up as an invalid Edwards
| point after an ECM stage 1 (16K FFT), i.e. a multithreaded FFT multiply
| occasionally returns a corrupted product.
|
| Method: run the *same deterministic sequence* of gwnum operations twice,
| once with num_threads=1 (reference) and once with num_threads=N, and
| compare a digest of every result.  gwnum is deterministic, so any
| difference is MT corruption, and the first differing iteration plus the
| operation variant localises the bug.  Inputs are fresh random values in
| [0, 2^n-1) (modulus 2^n-1, same as mprime ECM/PRP on M<n>).
|
| Variants exercised (all use fresh, non-FFTed inputs, exactly like
| mprime's ECM stage 1 does):
|   0 gwmul3      a*b            GWMUL_FFT_S12
|   1 gwmul3      a*b            GWMUL_FFT_S12 | STARTNEXTFFT
|   2 gwmulmuladd5 a*b+c*d       GWMUL_FFT_S1234
|   3 gwmulmuladd5 a*b+c*d       GWMUL_FFT_S1234 | STARTNEXTFFT
|   4 gwsubmul4   (a-b)*c        GWMUL_FFT_S123 | STARTNEXTFFT
|   5 chained:    t=a*b (STARTNEXTFFT); dst=t*c   (FFT_S2|STARTNEXTFFT)
|                  - the FFTed-source path that ECM's ed_extend/ed_add use
|
| Public gwnum API only (gwnum.h / giants.h).
|
| Build x86_64-windows:
|   g++ -O2 -o mtstress.exe mtstress.c -I<repo>/gwnum \
|       <work>/gwnum/lib-x86/gwnum.a <work>/gwnum/lib-x86/polymult.a \
|       -lstdc++ -lm -lpthread
| Build aarch64 static:
|   <zig> c++ -target aarch64-linux-musl -O2 -static -o mtstress-arm \
|       mtstress.c -I<repo>/gwnum <repo>/gwnum/gwnum.a \
|       <repo>/gwnum/polymult.a -lm -lpthread
|
| Usage: mtstress [--n EXP] [--threads N] [--iters K] [--seed S]
|                 [--op-list "0,1,2,3,4,5"] [--first-mismatch-only] [--quiet]
| Exit code: 0 = all digests matched, 1 = mismatch found, 2 = setup error.
+---------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "gwnum.h"
#include "cpuid.h"

static unsigned long	opt_n = 264169;		/* modulus 2^n-1 */
static int		opt_threads = 8;
static long		opt_iters = 20000;
static uint64_t		opt_seed = 0x123456789abcdef0ULL;
static int		opt_quiet = 0;
static int		opt_stop_first = 1;
static char		opt_ops[64] = "012345";

/* ---------------- deterministic LCG ---------------- */

static uint64_t lcg_state;
static uint32_t lcg (void)
{
	lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
	return ((uint32_t) (lcg_state >> 33));
}

/* ---------------- FNV-1a digest of a gwnum ---------------- */

static uint64_t fnv (uint64_t h, const void *data, size_t len)
{
	const unsigned char *p = (const unsigned char *) data;
	size_t i;
	for (i = 0; i < len; i++) {
		h ^= p[i];
		h *= 1099511628211ULL;
	}
	return h;
}

static int digest_err;		/* set if gwtogiant failed */

static uint64_t digest_gwnum (gwhandle *gw, gwnum x, giant buf)
{
	uint64_t h = 1469598103934665603ULL;
	int err = gwtogiant (gw, x, buf);
	if (err < 0) {
		digest_err = 1;
		return (h);
	}
	if (buf->sign == 0) return (fnv (h, "z", 1));
	if (buf->sign < 0) h = fnv (h, "-", 1);
	{
		int words = buf->sign < 0 ? -buf->sign : buf->sign;
		h = fnv (h, &words, sizeof (words));
	}
	h = fnv (h, buf->n, (size_t) (buf->sign < 0 ? -buf->sign : buf->sign) * sizeof (uint32_t));
	return (h);
}

/* ---------------- test context ---------------- */

static int	n_words;
static uint32_t	topmask;

typedef struct {
	gwhandle	gw;
	gwnum		a, b, c, d, t, dst;
	giant		gbuf;
} mtctx;

static void rand_giant (giant g, int words)
{
	int i;
	g->sign = words;
	for (i = 0; i < words; i++) g->n[i] = lcg ();
	g->n[words - 1] &= topmask;
	if (g->n[words - 1] == 0) g->n[words - 1] = 1;
}

/* Fill a gwnum with a fresh random value (NOT FFTed). */
static void rand_gwnum (mtctx *ctx, gwnum x)
{
	rand_giant (ctx->gbuf, n_words);
	gianttogw (&ctx->gw, ctx->gbuf, x);
}

/* ---------------- one pass ---------------- */

typedef struct {
	int		iter;
	int		variant;
	uint64_t	got;
	uint64_t	want;
} mismatch;

static int op_index (int iter, int nops)
{
	int which = iter % nops;
	return (opt_ops[which] - '0');
}

/* Run `iters` iterations at the given thread count.  mode 0 = fill `ref`
   with digests (reference pass); mode 1 = compare against `ref` and record
   mismatches.  Returns 0 on success, -1 on setup error. */
static int run_pass (mtctx *ctx, int threads, long iters, uint64_t *ref, int mode,
		     mismatch *first, int *num_mismatch)
{
	long i;
	int nops = (int) strlen (opt_ops);
	clock_t t0;

	gwinit (&ctx->gw);
	gwset_num_threads (&ctx->gw, threads);
	gwset_maxmulbyconst (&ctx->gw, 2);	/* GWMUL_MULBYCONST */
	/* Deterministic FFT choice: no benchmarks, single worker/one core would
	   change the jmptab selection between the two passes if left to default. */
	gwset_bench_cores (&ctx->gw, 1);
	gwset_bench_workers (&ctx->gw, 1);
	gwclear_use_benchmarks (&ctx->gw);
	gwset_safety_margin (&ctx->gw, 0.0);
	gwset_larger_fftlen_count (&ctx->gw, 0);
	gwset_minimum_fftlen (&ctx->gw, 0);
	gwset_num_threads (&ctx->gw, threads);
	if (gwsetup (&ctx->gw, 1.0, 2.0, (double) opt_n, -1.0)) {
		printf ("gwsetup failed\n");
		return (-1);
	}
	if (!opt_quiet) {
		char desc[256];
		gwfft_description (&ctx->gw, desc);
		printf ("  threads=%d  %s\n", threads, desc);
		fflush (stdout);
	}
	ctx->a = gwalloc (&ctx->gw);
	ctx->b = gwalloc (&ctx->gw);
	ctx->c = gwalloc (&ctx->gw);
	ctx->d = gwalloc (&ctx->gw);
	ctx->t = gwalloc (&ctx->gw);
	ctx->dst = gwalloc (&ctx->gw);
	ctx->gbuf = popg (&ctx->gw.gdata, n_words + 4);
	if (!ctx->a || !ctx->b || !ctx->c || !ctx->d || !ctx->t || !ctx->dst || !ctx->gbuf) {
		printf ("gwalloc failed\n");
		return (-1);
	}

	lcg_state = opt_seed;		/* identical input sequence for both passes */
	digest_err = 0;
	t0 = clock ();
	for (i = 0; i < iters; i++) {
		int v = op_index ((int) i, nops);
		uint64_t h = 0;
		int failed = 0;

		switch (v) {
		case 0:		/* a*b, no next-fft */
			rand_gwnum (ctx, ctx->a);
			rand_gwnum (ctx, ctx->b);
			gwmul3 (&ctx->gw, ctx->a, ctx->b, ctx->dst, GWMUL_FFT_S12);
			break;
		case 1:		/* a*b, start next fft */
			rand_gwnum (ctx, ctx->a);
			rand_gwnum (ctx, ctx->b);
			gwmul3 (&ctx->gw, ctx->a, ctx->b, ctx->dst,
				GWMUL_FFT_S12 | GWMUL_STARTNEXTFFT);
			break;
		case 2:		/* a*b + c*d, no next-fft */
			rand_gwnum (ctx, ctx->a);
			rand_gwnum (ctx, ctx->b);
			rand_gwnum (ctx, ctx->c);
			rand_gwnum (ctx, ctx->d);
			gwmulmuladd5 (&ctx->gw, ctx->a, ctx->b, ctx->c, ctx->d,
				      ctx->dst, GWMUL_FFT_S1234);
			break;
		case 3:		/* a*b + c*d, start next fft */
			rand_gwnum (ctx, ctx->a);
			rand_gwnum (ctx, ctx->b);
			rand_gwnum (ctx, ctx->c);
			rand_gwnum (ctx, ctx->d);
			gwmulmuladd5 (&ctx->gw, ctx->a, ctx->b, ctx->c, ctx->d,
				      ctx->dst, GWMUL_FFT_S1234 | GWMUL_STARTNEXTFFT);
			break;
		case 4:		/* (a-b)*c, start next fft */
			rand_gwnum (ctx, ctx->a);
			rand_gwnum (ctx, ctx->b);
			rand_gwnum (ctx, ctx->c);
			gwsubmul4 (&ctx->gw, ctx->a, ctx->b, ctx->c, ctx->dst,
				   GWMUL_FFT_S123 | GWMUL_STARTNEXTFFT);
			break;
		case 5:		/* chained: t = a*b (ffted); dst = t*c */
			rand_gwnum (ctx, ctx->a);
			rand_gwnum (ctx, ctx->b);
			rand_gwnum (ctx, ctx->c);
			gwmul3 (&ctx->gw, ctx->a, ctx->b, ctx->t,
				GWMUL_FFT_S12 | GWMUL_STARTNEXTFFT);
			gwmul3 (&ctx->gw, ctx->t, ctx->c, ctx->dst,
				GWMUL_FFT_S2 | GWMUL_STARTNEXTFFT);
			break;
		default:
			printf ("bad variant %d\n", v);
			return (-1);
		}

		if (gw_test_for_error (&ctx->gw)) failed = 1;
		h = digest_gwnum (&ctx->gw, ctx->dst, ctx->gbuf);
		if (digest_err) failed = 1;

		if (mode == 1) {
			if (h != ref[i] || failed) {
				(*num_mismatch)++;
				if (*num_mismatch <= 20) {
					printf ("  MISMATCH iter %ld variant %d op=%s err=%d got=%016llx want=%016llx\n",
						i, v, "?", failed,
						(unsigned long long) h,
						(unsigned long long) ref[i]);
					fflush (stdout);
				}
				if (*num_mismatch == 1 && first) {
					first->iter = (int) i;
					first->variant = v;
					first->got = h;
					first->want = ref[i];
				}
				if (opt_stop_first) break;
			}
		} else {
			ref[i] = h | (failed ? 0x1000000000000000ULL : 0);
		}
	}
	if (!opt_quiet) {
		double secs = (double) (clock () - t0) / (double) CLOCKS_PER_SEC;
		printf ("  threads=%d done %ld iters in %.2f s (%.1f iters/s)%s\n",
			threads, iters, secs, secs > 0 ? (double) iters / secs : 0.0,
			digest_err ? "  [gwtogiant errors!]" : "");
		fflush (stdout);
	}

	gwfree (&ctx->gw, ctx->a);
	gwfree (&ctx->gw, ctx->b);
	gwfree (&ctx->gw, ctx->c);
	gwfree (&ctx->gw, ctx->d);
	gwfree (&ctx->gw, ctx->t);
	gwfree (&ctx->gw, ctx->dst);
	pushg (&ctx->gw.gdata, 1);
	gwdone (&ctx->gw);
	memset (ctx, 0, sizeof (*ctx));
	return (0);
}

int main (int argc, char **argv)
{
	mtctx ctx;
	uint64_t *ref;
	long i;
	int arg, num_mismatch = 0;
	mismatch first;

	memset (&first, 0, sizeof (first));

	for (arg = 1; arg < argc; arg++) {
		if (!strcmp (argv[arg], "--n") && arg + 1 < argc) opt_n = strtoul (argv[++arg], NULL, 10);
		else if (!strcmp (argv[arg], "--threads") && arg + 1 < argc) opt_threads = atoi (argv[++arg]);
		else if (!strcmp (argv[arg], "--iters") && arg + 1 < argc) opt_iters = atol (argv[++arg]);
		else if (!strcmp (argv[arg], "--seed") && arg + 1 < argc) opt_seed = strtoull (argv[++arg], NULL, 0);
		else if (!strcmp (argv[arg], "--op-list") && arg + 1 < argc) { strncpy (opt_ops, argv[++arg], sizeof (opt_ops) - 1); opt_ops[sizeof (opt_ops) - 1] = 0; }
		else if (!strcmp (argv[arg], "--first-mismatch-only")) opt_stop_first = 1;
		else if (!strcmp (argv[arg], "--all-mismatches")) opt_stop_first = 0;
		else if (!strcmp (argv[arg], "--quiet")) opt_quiet = 1;
		else if (!strcmp (argv[arg], "--help")) {
			printf ("usage: mtstress [--n EXP] [--threads N] [--iters K] [--seed S]\n"
				"                [--op-list 012345] [--first-mismatch-only|--all-mismatches] [--quiet]\n");
			return (0);
		} else {
			printf ("unknown option %s\n", argv[arg]);
			return (2);
		}
	}

	if (opt_iters <= 0) { printf ("bad iters\n"); return (2); }
	n_words = (int) (opt_n >> 5) + 1;
	topmask = (opt_n & 31) ? ((1u << (opt_n & 31)) - 1) : 0xffffffffu;

	printf ("mtstress: n=%lu (2^n-1), threads=%d, iters=%ld, seed=%016llx, ops=%s\n",
		opt_n, opt_threads, opt_iters, (unsigned long long) opt_seed, opt_ops);
	guessCpuType ();
	/* Make the x86 host build behave exactly like the ARM64 build: the
	   converted kernels only implement the 128-bit SSE2 family, so the
	   implementation/FFT selection must see the same CPU_FLAGS as ARM64
	   (SSE..SSE42 only).  Without this the x86 build selects a different
	   FFT layout and we would not be testing the phone's code path. */
	CPU_FLAGS &= ~(CPU_AVX | CPU_FMA3 | CPU_FMA4 | CPU_AVX2 |
		       CPU_AVX512F | CPU_AVX512PF | CPU_AVX512DQ | CPU_AVX512VL |
		       CPU_PREFETCHW | CPU_PREFETCHWT1);
	/* ARM64 also reports architecture "other" with zero family/model/signature
	   (gwnum's jmptable selection looks at the architecture). */
	CPU_ARCHITECTURE = CPU_ARCHITECTURE_OTHER;
	CPU_FAMILY = 0;
	CPU_MODEL = 0;
	CPU_SIGNATURE = 0;
	printf ("CPU: %s, %.0f MHz, %u cores, %u hyperthreads, flags=0x%08X (ARM64-equivalent)\n",
		CPU_BRAND, CPU_SPEED, CPU_CORES, CPU_HYPERTHREADS, CPU_FLAGS);
	printf ("phase 1: reference, threads=1\n");
	fflush (stdout);

	ref = (uint64_t *) malloc ((size_t) opt_iters * sizeof (uint64_t));
	if (!ref) { printf ("out of memory\n"); return (2); }

	memset (&ctx, 0, sizeof (ctx));
	if (run_pass (&ctx, 1, opt_iters, ref, 0, NULL, &num_mismatch) != 0) return (2);

	if (opt_threads > 1) {
		printf ("phase 2: threads=%d\n", opt_threads);
		fflush (stdout);
		num_mismatch = 0;
		if (run_pass (&ctx, opt_threads, opt_iters, ref, 1, &first, &num_mismatch) != 0) return (2);
	}

	if (num_mismatch == 0) {
		printf ("RESULT: PASS - %ld iterations bit-identical at %d thread(s)\n",
			opt_iters, opt_threads);
		free (ref);
		return (0);
	}
	printf ("RESULT: CORRUPTION - %d mismatching iteration(s); first at iter %d (variant %d) got=%016llx want=%016llx\n",
		num_mismatch, first.iter, first.variant,
		(unsigned long long) first.got, (unsigned long long) first.want);
	free (ref);
	return (1);
}
