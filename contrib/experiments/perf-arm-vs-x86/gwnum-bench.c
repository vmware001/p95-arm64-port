/* gwnum FFT throughput microbenchmark (single thread, fixed FFT size).
 *
 * Measures the cost of one gwnum multiply/square at a given exponent using only
 * the public gwnum API, so the identical C source can be built against:
 *   - our asm->C converted gwnum for x86_64 and aarch64 (same kernels)
 *   - the official x86_64 gwnum (hand written asm) for comparison
 *
 * The x86 build forces the CPU-flags/architecture profile of the ARM64 build, so
 * both builds select the same 128-bit FFT layout (e.g. "type-1 FFT length 16K,
 * Pass1=64, Pass2=256, clm=4").
 *
 * Usage: gwnum-bench [p] [iterations] [mul|square]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "gwnum.h"
#include "cpuid.h"

#if defined(_WIN32)
#include <windows.h>
static double now_sec(void)
{
	static LARGE_INTEGER freq;
	LARGE_INTEGER c;
	if (freq.QuadPart == 0) QueryPerformanceFrequency(&freq);
	QueryPerformanceCounter(&c);
	return (double) c.QuadPart / (double) freq.QuadPart;
}
static uint64_t read_cycles(void)
{
	/* QueryThreadCycleTime returns the real CPU cycles this thread executed,
	   so it also works with turbo/boost (unlike the fixed-frequency TSC). */
	ULONG64 c = 0;
	if (!QueryThreadCycleTime(GetCurrentThread(), &c)) return __rdtsc();
	return (uint64_t) c;
}
#elif defined(__x86_64__)
#include <x86intrin.h>
static double now_sec(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}
static uint64_t read_cycles(void) { return __rdtsc(); }
#else
static double now_sec(void)
{
	struct timespec ts;
	clock_gettime(CLOCK_MONOTONIC, &ts);
	return ts.tv_sec + ts.tv_nsec * 1e-9;
}
static uint64_t read_cycles(void) { return 0; }
#endif

static uint32_t rng_state = 123456789u;
static uint32_t rng(void)
{
	rng_state = rng_state * 1664525u + 1013904223u;
	return rng_state;
}

int main(int argc, char **argv)
{
	unsigned long p = 264169;
	long iters = 2000;
	int mode = 0;			/* 0 = multiply a*b, 1 = square */
	int force_128bit = 1;		/* x86: force the ARM64-equivalent CPU profile */
	gwhandle gw;
	gwnum a, b, c;
	giant g;
	int words, i, err;
	double t0, t1, secs;
	uint64_t c0 = 0, c1 = 0;
	char desc[256];

	if (argc > 1) p = strtoul(argv[1], NULL, 10);
	if (argc > 2) iters = atol(argv[2]);
	if (argc > 3) mode = !strcmp(argv[3], "square");
	if (argc > 4 && !strcmp(argv[4], "native")) force_128bit = 0;

#if defined(__x86_64__) || defined(_M_X64)
	/* Make the x86 build select the same 128-bit code paths and FFT layout as
	   the ARM64 build (see contrib/experiments/mt-fix/mtstress.c).  With
	   "native" as the 4th argument, use the CPU's native ISA instead. */
	guessCpuType();
	if (force_128bit) {
		CPU_FLAGS &= ~(CPU_AVX | CPU_FMA3 | CPU_FMA4 | CPU_AVX2 |
			       CPU_AVX512F | CPU_AVX512PF | CPU_AVX512DQ | CPU_AVX512VL |
			       CPU_PREFETCHW | CPU_PREFETCHWT1);
		CPU_ARCHITECTURE = CPU_ARCHITECTURE_OTHER;
		CPU_FAMILY = CPU_MODEL = 0;
		CPU_SIGNATURE = 0;
	}
	printf("cpu_flags=0x%08X (%s)\n", CPU_FLAGS, force_128bit ? "forced 128-bit" : "native");
#endif

	gwinit(&gw);
	gwset_num_threads(&gw, 1);
	gwset_bench_cores(&gw, 1);
	gwset_bench_workers(&gw, 1);
	gwclear_use_benchmarks(&gw);
	gwset_safety_margin(&gw, 0.0);
	gwset_larger_fftlen_count(&gw, 0);
	gwset_minimum_fftlen(&gw, 0);
	err = gwsetup(&gw, 1.0, 2, p, -1);
	if (err) { printf("gwsetup failed: %d\n", err); return 2; }
	gwfft_description(&gw, desc);

	a = gwalloc(&gw);
	b = gwalloc(&gw);
	c = gwalloc(&gw);
	if (!a || !b || !c) { printf("gwalloc failed\n"); return 2; }

	words = (int) (p >> 5) + 2;
	g = popg(&gw.gdata, words);
	g->sign = words;
	for (i = 0; i < words; i++) g->n[i] = rng();
	gianttogw(&gw, g, a);
	for (i = 0; i < words; i++) g->n[i] = rng();
	gianttogw(&gw, g, b);
	pushg(&gw.gdata, 1);

	/* Warm-up */
	for (i = 0; i < (iters > 20 ? 20 : iters); i++) {
		if (mode) gwsquare2(&gw, a, c, 0);
		else gwmul3(&gw, a, b, c, 0);
	}

	c0 = read_cycles();
	t0 = now_sec();
	for (i = 0; i < iters; i++) {
		if (mode) gwsquare2(&gw, a, c, 0);
		else gwmul3(&gw, a, b, c, 0);
	}
	t1 = now_sec();
	c1 = read_cycles();
	secs = t1 - t0;

	printf("p=%lu iters=%ld mode=%s fft=\"%s\"\n", p, iters, mode ? "square" : "mul", desc);
	printf("  time    = %.6f s total, %.3f us/op, %.1f ops/s\n",
	       secs, secs * 1e6 / iters, iters / secs);
	if (c1 > c0) {
		double tsc_hz = (double) (c1 - c0) / secs;
		printf("  cycles  = %.0f cycles/op (tsc %.2f MHz)\n",
		       (double) (c1 - c0) / iters, tsc_hz / 1e6);
	} else {
#if defined(__aarch64__)
		/* AArch64 user space has no rdtsc; the caller supplies the CPU clock
		   (RK3528 runs at a fixed 2016 MHz with the performance governor). */
		printf("  cycles  = %.0f cycles/op (assuming %.0f MHz CPU clock)\n",
		       secs * 1e6 / iters * 2016.0, 2016.0);
#else
		printf("  cycles  = n/a\n");
#endif
	}

	gwfree(&gw, a);
	gwfree(&gw, b);
	gwfree(&gw, c);
	gwdone(&gw);
	return 0;
}
