// fftbench_asm.c -- gwnum squaring benchmark + bit-exactness digest.
// Usage: fftbench_asm <exponent> <iters>
// Output (stdout, stable format, one field per line):
//   fftlen=<n>
//   timeA_ms=<x>            (timed loop, operand = 3.0)
//   timeB_ms=<x>            (timed loop, operand = full-size pseudo-random)
//   digest=<16 hex digits>  (FNV-1a-64 over raw FFT doubles + giant limbs of
//                            several full-size multiply results)
// Env: PROFILE_US=<n> -> additionally sample PCs during an extra timed loop
//      and dump them to stderr.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include "gwnum.h"

void prof_start(int usec);
void prof_stop(void);
void prof_dump(const char *label);
void prof_reset(void);

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static uint64_t fnv_buf(uint64_t h, const void *p, size_t nbytes) {
    const unsigned char *b = (const unsigned char *) p;
    for (size_t i = 0; i < nbytes; i++) { h ^= b[i]; h *= 1099511628211ULL; }
    return h;
}

static uint64_t rng_state;
static uint32_t rng_next(void) {
    rng_state = rng_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return (uint32_t) (rng_state >> 32);
}

int main(int argc, char **argv) {
    if (argc < 3) { fprintf(stderr, "usage: fftbench_asm <exp> <iters>\n"); return 2; }
    unsigned long n = strtoul(argv[1], 0, 10);
    long iters = strtol(argv[2], 0, 10);

    gwhandle gw;
    gwinit(&gw);
    if (gwsetup(&gw, 1, 2, n, -1)) { fprintf(stderr, "gwsetup failed\n"); return 1; }
    unsigned long fftlen = gwfftlen(&gw);
    printf("fftlen=%lu\n", fftlen);

    gwnum x = gwalloc(&gw), z = gwalloc(&gw), y = gwalloc(&gw);
    dbltogw(&gw, 3.0, x);

    // ---- timing A: constant small operand (matches historical baseline) ----
    for (int i = 0; i < 50; i++) gwmul3(&gw, x, x, z, 0);
    double t0 = now_sec();
    for (long i = 0; i < iters; i++) gwmul3(&gw, x, x, z, 0);
    double dtA = now_sec() - t0;
    printf("timeA_ms=%.4f\n", 1000.0 * dtA / iters);

    // ---- build a full-size pseudo-random operand, time that too ----
    int nwords = (int) (n / 32);
    uint32_t topmask = (uint32_t) ((n % 32) ? ((1ULL << (n % 32)) - 1) : 0xffffffffu);
    giant gx = allocgiant(nwords + 12);
    giant gz = allocgiant(nwords + 12);
    giant gy = allocgiant(nwords + 12);

    uint64_t digest = 1469598103934665603ULL;

    rng_state = 0x123456789abcdefULL;
    for (int k = 0; k < 6; k++) {
        for (int j = 0; j < nwords; j++) gx->n[j] = rng_next();
        if (n % 32) gx->n[nwords] = rng_next() & topmask;
        else gx->n[nwords - 1] &= 0xffffffffu;
        int nsig = (n % 32) ? nwords + 1 : nwords;
        while (nsig > 1 && gx->n[nsig - 1] == 0) nsig--;
        gx->sign = nsig;

        gianttogw(&gw, gx, x);
        if (k == 0) {
            for (int i = 0; i < 20; i++) gwmul3(&gw, x, x, z, 0);
            double tb0 = now_sec();
            for (long i = 0; i < iters; i++) gwmul3(&gw, x, x, z, 0);
            double dtB = now_sec() - tb0;
            printf("timeB_ms=%.4f\n", 1000.0 * dtB / iters);
        }
        gwmul3(&gw, x, x, z, 0);   // z = x^2 mod 2^n-1
        digest = fnv_buf(digest, (const void *) z, sizeof(double) * fftlen);
        gwtogiant(&gw, z, gz);
        digest = fnv_buf(digest, gz->n, sizeof(uint32_t) * (size_t) (gz->sign < 0 ? -gz->sign : gz->sign));
        // cross product x * z
        gwmul3(&gw, x, z, y, 0);
        digest = fnv_buf(digest, (const void *) y, sizeof(double) * fftlen);
        gwtogiant(&gw, y, gy);
        digest = fnv_buf(digest, gy->n, sizeof(uint32_t) * (size_t) (gy->sign < 0 ? -gy->sign : gy->sign));
    }
    printf("digest=%016llx\n", (unsigned long long) digest);

    // ---- optional PC sampling of the hot loop ----
    {
        const char *pus = getenv("PROFILE_US");
        if (pus && atoi(pus) > 0) {
            gianttogw(&gw, gx, x);   // full-size operand
            prof_reset();
            prof_start(atoi(pus));
            for (long i = 0; i < iters; i++) gwmul3(&gw, x, x, z, 0);
            prof_stop();
            prof_dump("bench");
            printf("profile_digest=%016llx\n",
                   (unsigned long long) fnv_buf(1469598103934665603ULL, (const void *) z, sizeof(double) * fftlen));
        }
    }

    freegiant(gx); freegiant(gz); freegiant(gy);
    gwdone(&gw);
    return 0;
}
