// copytest.c -- differential bit-exactness + micro-benchmark for the
// hand written NEON pass1 copy kernel.
// Usage: copytest [reps]
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <time.h>

void r4_pass1_copy_asm(void *src, void *dst, long src_mid, long dst_mid,
                       long src_z, long dst_z, long mid_cnt, long z_cnt);

#define SZ (1 << 20)

static double now_sec(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

// reference: exactly the semantics documented for r4_pass1_copy_asm
static void __attribute__((noinline)) ref_copy(unsigned char *src, unsigned char *dst,
                     long src_mid, long dst_mid, long src_z, long dst_z,
                     long mid_cnt, long z_cnt) {
    for (long z = 0; z < z_cnt; z++) {
        for (long y = 0; y < mid_cnt; y++) {
            memcpy(dst, src, 512);
            src += src_mid;
            dst += dst_mid;
        }
        src += src_z;
        dst += dst_z;
    }
}

int main(int argc, char **argv) {
    long reps = (argc > 1) ? strtol(argv[1], 0, 10) : 20000;
    unsigned char *src = aligned_alloc(64, SZ);
    unsigned char *d1 = aligned_alloc(64, SZ);
    unsigned char *d2 = aligned_alloc(64, SZ);
    if (!src || !d1 || !d2) { fprintf(stderr, "alloc failed\n"); return 1; }
    uint64_t s = 12345;
    for (int i = 0; i < SZ; i++) {
        s = s * 6364136223846793005ULL + 1442695040888963407ULL;
        src[i] = (unsigned char) (s >> 33);
    }

    // geometry used by the 16K FFT pass 1 kernel (clm = 4):
    //   scratch -> fft : src_mid = clm*128 = 512, dst_mid = blkdst = 8320,
    //                    src_z = 128, dst_z = 0, mid = 8, z = 2 (cachelines/8)
    //   fft -> scratch : swapped
    struct { long sm, dm, sz, dz, mid, z; } cases[2] = {
        {512, 8320, 128, 0, 8, 2},
        {8320, 512, 0, 128, 8, 2},
    };
    int bad = 0;
    for (int c = 0; c < 2; c++) {
        memset(d1, 0, SZ);
        memset(d2, 0, SZ);
        r4_pass1_copy_asm(src, d1, cases[c].sm, cases[c].dm, cases[c].sz,
                          cases[c].dz, cases[c].mid, cases[c].z);
        ref_copy(src, d2, cases[c].sm, cases[c].dm, cases[c].sz, cases[c].dz,
                 cases[c].mid, cases[c].z);
        if (memcmp(d1, d2, SZ) != 0) {
            bad++;
            for (int i = 0; i < SZ; i++)
                if (d1[i] != d2[i]) {
                    printf("case %d MISMATCH at byte %d: asm=%02x ref=%02x\n",
                           c, i, d1[i], d2[i]);
                    break;
                }
        } else {
            printf("case %d bit-exact ok\n", c);
        }
    }
    if (bad) { printf("COPYTEST FAIL\n"); return 1; }

    // micro benchmark: many repeats of the exact kernel geometry
    const int iters = (int) reps;
    volatile long sink = 0;
    struct timespec t0, t1;
    long long na, nr;
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < iters; i++) {
        r4_pass1_copy_asm(src, d1, 512, 8320, 128, 0, 8, 2);
        sink += d1[i & 1023];
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    na = (long long) (t1.tv_sec - t0.tv_sec) * 1000000000LL + (t1.tv_nsec - t0.tv_nsec);
    clock_gettime(CLOCK_MONOTONIC, &t0);
    for (int i = 0; i < iters; i++) {
        ref_copy(src, d2, 512, 8320, 128, 0, 8, 2);
        sink += d2[i & 1023];
    }
    clock_gettime(CLOCK_MONOTONIC, &t1);
    nr = (long long) (t1.tv_sec - t0.tv_sec) * 1000000000LL + (t1.tv_nsec - t0.tv_nsec);
    printf("asm  %lld ns/call\n", na / iters);
    printf("ref  %lld ns/call\n", nr / iters);
    printf("speedup %.3f x  (sink=%ld)\n", (double) nr / (double) na, sink);
    printf("bitsum %llx\n",
           (unsigned long long) ({ uint64_t h = 0; for (int i = 0; i < SZ; i++) h = h * 131 + d1[i]; h; }));
    printf("COPYTEST PASS\n");
    return 0;
}
