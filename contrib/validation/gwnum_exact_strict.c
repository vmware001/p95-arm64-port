/* gwnum_exact_strict.c -- bit-exact validation of gwnum arithmetic against GMP.
 *
 * Purpose: after any optimisation of the gwnum FFT kernels (compiler flags,
 * assembly, layout changes) we must prove that the arithmetic is still EXACT.
 * gwnum's FFT multiply is exact by design (the FFT size is chosen so that the
 * roundoff is far below the carry threshold), therefore for every input pair
 * the recovered integer must equal the true value a*b mod 2^n-1 exactly.
 * Any difference - even in the last bit - means the optimisation broke
 * correctness and would silently corrupt ECM/PRP results.
 *
 * Three independent checks per iteration:
 *   mul  : gwmul3(x, y)     == (a*b)      mod 2^n-1
 *   add  : gwadd3o(x, y)    == (a+b)      mod 2^n-1
 *   sqr  : gwmul3(x, x)     == (a*a)      mod 2^n-1
 * plus deliberately correlated "worst case" inputs (a = b = 2^n-1-1, a = b = 1),
 * because that is where the FFT roundoff is largest.
 *
 * Small exponents are used (GMP is cheap there) plus one large "production"
 * size so that the same FFT kernels that the farm actually runs are exercised.
 *
 * Build:
 *   zig cc -target aarch64-linux-musl -O2 -I<gwnum dir> -I<deps>/include \
 *          -c gwnum_exact_strict.c -o gwnum_exact_strict.o
 *   zig cc -target aarch64-linux-musl -static -no-pie -s -o gwnum_exact_strict \
 *          gwnum_exact_strict.o libgwnum.a libpolymult.a libgmp.a -lm -lpthread
 *
 * Usage: gwnum_exact_strict <exponent> <iterations>
 * Exit code 0 = all exact, 3 = mismatch found.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

#include "gwnum.h"
#include "giants.h"
#include <gmp.h>

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static uint32_t *g_limbs(giant g) { return g->n; }
static long g_nlimbs(giant g) { return (long) (g->sign < 0 ? -g->sign : g->sign); }

/* mpz -> giant (little-endian 32-bit limbs), zero-filling the whole buffer */
static void mpz_to_giant(giant g, unsigned long cap, mpz_srcptr z)
{
    memset(g_limbs(g), 0, (size_t) cap * 4);
    size_t count = 0;
    mpz_export(g_limbs(g), &count, -1, sizeof(uint32_t), 0, 0, z);
    g->sign = (int) count;
}

/* compare gwnum (converted through its giant) with the exact mpz value */
static int cmp_gw_mpz(gwhandle *gw, gwnum x, giant g, unsigned long cap, mpz_srcptr exact,
                      const char *what, unsigned long exp, long iter)
{
    int rc = gwtogiant(gw, x, g);
    if (rc) {
        fprintf(stderr, "FAIL %s exp=%lu iter=%ld: gwtogiant rc=%d\n", what, exp, iter, rc);
        return 1;
    }
    long n = g_nlimbs(g);
    mpz_t got;
    mpz_init(got);
    mpz_import(got, (size_t) n, -1, sizeof(uint32_t), 0, 0, g_limbs(g));
    int bad = (mpz_cmp(got, exact) != 0);
    if (bad) {
        fprintf(stderr, "FAIL %s exp=%lu iter=%ld: got %ld limbs vs exact %ld limbs\n",
                what, exp, iter, n, (long) mpz_size(exact));
    }
    mpz_clear(got);
    return bad;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "usage: gwnum_exact_strict <exp> <iters>\n");
        return 2;
    }
    unsigned long n = strtoul(argv[1], 0, 10);
    long iters = strtol(argv[2], 0, 10);
    unsigned long cap = n / 32 + 64;

    gwhandle gw;
    gwinit(&gw);
    if (gwsetup(&gw, 1, 2, n, -1)) { fprintf(stderr, "gwsetup failed\n"); return 1; }

    gwnum x = gwalloc(&gw), y = gwalloc(&gw), z = gwalloc(&gw);
    giant gx = allocgiant(cap), gy = allocgiant(cap), gz = allocgiant(cap);
    if (!gx || !gy || !gz) { fprintf(stderr, "allocgiant failed\n"); return 1; }

    mpz_t M, X, Y, T, S;
    mpz_init(M); mpz_set_ui(M, 1); mpz_mul_2exp(M, M, n); mpz_sub_ui(M, M, 1);
    mpz_inits(X, Y, T, S, NULL);

    /* deterministic PRNG so runs are reproducible */
    gmp_randstate_t rs;
    gmp_randinit_mt(rs);
    gmp_randseed_ui(rs, 0x5EED1234);

    long fails = 0, checks = 0;
    double t0 = now_sec();

    for (long i = 0; i < iters; i++) {
        int mode = (int) (i % 4);
        if (mode == 0) {                     /* random a,b */
            mpz_urandomb(X, rs, n); mpz_urandomb(Y, rs, n);
        } else if (mode == 1) {              /* correlated worst case */
            mpz_set(X, M); mpz_sub_ui(X, X, 1); mpz_set(Y, X);
        } else if (mode == 2) {              /* a = b = 1 */
            mpz_set_ui(X, 1); mpz_set_ui(Y, 1);
        } else {                             /* random near modulus */
            mpz_urandomb(X, rs, n); mpz_set_ui(Y, 1); mpz_mul_2exp(Y, Y, n / 2);
            mpz_add(Y, Y, M); mpz_sub_ui(Y, Y, 1);
        }

        mpz_to_giant(gx, cap, X);
        mpz_to_giant(gy, cap, Y);
        gianttogw(&gw, gx, x);
        gianttogw(&gw, gy, y);

        /* --- multiply: z = x*y mod M --- */
        gwmul3(&gw, x, y, z, 0);
        mpz_mul(T, X, Y); mpz_mod(T, T, M);
        fails += cmp_gw_mpz(&gw, z, gz, cap, T, "mul", n, i); checks++;

        /* --- add: z = x+y mod M --- */
        mpz_to_giant(gx, cap, X);
        mpz_to_giant(gy, cap, Y);
        gianttogw(&gw, gx, x);
        gianttogw(&gw, gy, y);
        gwadd3o(&gw, x, y, z, 0);
        mpz_add(T, X, Y); mpz_mod(T, T, M);
        fails += cmp_gw_mpz(&gw, z, gz, cap, T, "add", n, i); checks++;

        /* --- square: z = x^2 mod M --- */
        mpz_to_giant(gx, cap, X);
        gianttogw(&gw, gx, x);
        gwmul3(&gw, x, x, z, 0);
        mpz_mul(T, X, X); mpz_mod(T, T, M);
        fails += cmp_gw_mpz(&gw, z, gz, cap, T, "sqr", n, i); checks++;

        /* --- fully-FFTed multiply (different code path: both args FFTed) --- */
        mpz_to_giant(gx, cap, X);
        mpz_to_giant(gy, cap, Y);
        gianttogw(&gw, gx, x);
        gianttogw(&gw, gy, y);
        gwfft(&gw, x, x);
        gwfft(&gw, y, y);
        gwmul3(&gw, x, y, z, GWMUL_PRESERVE_S1 | GWMUL_PRESERVE_S2);
        mpz_mul(T, X, Y); mpz_mod(T, T, M);
        fails += cmp_gw_mpz(&gw, z, gz, cap, T, "fftmul", n, i); checks++;

        if (fails) break;
    }
    double dt = now_sec() - t0;

    printf("RESULT exp=%lu fftlen=%lu iters=%ld checks=%ld fails=%ld time=%.2fs %s\n",
           n, gwfftlen(&gw), iters, checks, fails, dt, fails ? "MISMATCH" : "ALL_EXACT");

    mpz_clears(M, X, Y, T, S, NULL);
    gmp_randclear(rs);
    freegiant(gx); freegiant(gy); freegiant(gz);
    gwdone(&gw);
    return fails ? 3 : 0;
}
