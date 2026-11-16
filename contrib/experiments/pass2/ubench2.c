// ubench2.c -- Cortex-A53 microbenchmark v2: dependency chains broken correctly.
// Measurements: DP FP latency/throughput, L1D load/store throughput, memory levels.
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

#define N 1000000
#define GHZ 2.016

// latency: one dependent chain of `insn` on v0
#define LAT(name, insn)                                                  \
static void name(double *out)                                            \
{                                                                        \
    double acc = 1.0, c = 1e-12;                                         \
    __asm__ volatile(                                                    \
        "fmov d16, %d[cin]\n" "fmov d0, %d[acc]\n"                       \
        "mov w9, %w[n]\n" "1:\n" insn " v0.2d, v0.2d, v16.2d\n"          \
        "subs w9, w9, #1\n" "b.ne 1b\n" "fmov %d[o], d0\n"               \
        : [o] "+w"(acc) : [acc] "w"(acc), [cin] "w"(c), [n] "r"(N)       \
        : "v0", "v16", "w9", "cc");                                      \
    out[0] = acc;                                                        \
}

// throughput: 8 independent chains
#define TP(name, insn)                                                   \
static void name(double *out)                                            \
{                                                                        \
    double acc = 1.0, c = 1e-12;                                         \
    __asm__ volatile(                                                    \
        "fmov d16, %d[cin]\n" "fmov d0, %d[acc]\n"                       \
        "fmov d1, d0\n" "fmov d2, d0\n" "fmov d3, d0\n"                  \
        "fmov d4, d0\n" "fmov d5, d0\n" "fmov d6, d0\n" "fmov d7, d0\n"  \
        "mov w9, %w[n]\n" "1:\n"                                         \
        insn " v0.2d, v0.2d, v16.2d\n" insn " v1.2d, v1.2d, v16.2d\n"    \
        insn " v2.2d, v2.2d, v16.2d\n" insn " v3.2d, v3.2d, v16.2d\n"    \
        insn " v4.2d, v4.2d, v16.2d\n" insn " v5.2d, v5.2d, v16.2d\n"    \
        insn " v6.2d, v6.2d, v16.2d\n" insn " v7.2d, v7.2d, v16.2d\n"    \
        "subs w9, w9, #1\n" "b.ne 1b\n"                                  \
        "fadd d0, d0, d1\n" "fadd d2, d2, d3\n" "fadd d4, d4, d5\n"      \
        "fadd d6, d6, d7\n" "fadd d0, d0, d2\n" "fadd d4, d4, d6\n"      \
        "fadd d0, d0, d4\n" "fmov %d[o], d0\n"                           \
        : [o] "+w"(acc) : [acc] "w"(acc), [cin] "w"(c), [n] "r"(N)       \
        : "v0","v1","v2","v3","v4","v5","v6","v7","v16","w9","cc");      \
    out[0] = acc;                                                        \
}

LAT(lat_fadd, "fadd")
LAT(lat_fmul, "fmul")
LAT(lat_fmla, "fmla")
TP(tp_fadd, "fadd")
TP(tp_fmul, "fmul")
TP(tp_fmla, "fmla")

// 16 loads per iter, 4 independent accumulators (chain length 4)
static void ld_q(double *out)
{
    static double buf[8192] __attribute__((aligned(64)));
    for (int i = 0; i < 8192; i++) buf[i] = i * 1e-9;
    double acc = 0;
    __asm__ volatile(
        "mov x9, %[b]\n" "mov w10, %w[n]\n"
        "movi v28.2d, #0\n" "movi v29.2d, #0\n" "movi v30.2d, #0\n" "movi v31.2d, #0\n"
        "1:\n"
        "ldr q0, [x9]\n" "ldr q1, [x9, #16]\n" "ldr q2, [x9, #32]\n" "ldr q3, [x9, #48]\n"
        "ldr q4, [x9, #64]\n" "ldr q5, [x9, #80]\n" "ldr q6, [x9, #96]\n" "ldr q7, [x9, #112]\n"
        "ldr q8, [x9, #128]\n" "ldr q9, [x9, #144]\n" "ldr q10, [x9, #160]\n" "ldr q11, [x9, #176]\n"
        "ldr q12, [x9, #192]\n" "ldr q13, [x9, #208]\n" "ldr q14, [x9, #224]\n" "ldr q15, [x9, #240]\n"
        "fadd v28.2d, v28.2d, v0.2d\n" "fadd v28.2d, v28.2d, v4.2d\n"
        "fadd v28.2d, v28.2d, v8.2d\n" "fadd v28.2d, v28.2d, v12.2d\n"
        "fadd v29.2d, v29.2d, v1.2d\n" "fadd v29.2d, v29.2d, v5.2d\n"
        "fadd v29.2d, v29.2d, v9.2d\n" "fadd v29.2d, v29.2d, v13.2d\n"
        "fadd v30.2d, v30.2d, v2.2d\n" "fadd v30.2d, v30.2d, v6.2d\n"
        "fadd v30.2d, v30.2d, v10.2d\n" "fadd v30.2d, v30.2d, v14.2d\n"
        "fadd v31.2d, v31.2d, v3.2d\n" "fadd v31.2d, v31.2d, v7.2d\n"
        "fadd v31.2d, v31.2d, v11.2d\n" "fadd v31.2d, v31.2d, v15.2d\n"
        "subs w10, w10, #1\n" "b.ne 1b\n"
        "fadd v28.2d, v28.2d, v29.2d\n" "fadd v30.2d, v30.2d, v31.2d\n"
        "fadd v28.2d, v28.2d, v30.2d\n" "fmov %d[o], d28\n"
        : [o] "+w"(acc) : [b] "r"(buf), [n] "r"(N)
        : "v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11",
          "v12","v13","v14","v15","v28","v29","v30","v31","x9","w10","cc");
    out[0] = acc;
}

// 16 stores per iter + 1 load
static void st_q(double *out)
{
    static double buf[8192] __attribute__((aligned(64)));
    double acc = 0;
    __asm__ volatile(
        "mov x9, %[b]\n" "mov w10, %w[n]\n"
        "movi v0.2d, #0\n" "movi v1.2d, #0\n" "movi v2.2d, #0\n" "movi v3.2d, #0\n"
        "1:\n"
        "str q0, [x9]\n" "str q1, [x9, #16]\n" "str q2, [x9, #32]\n" "str q3, [x9, #48]\n"
        "str q0, [x9, #64]\n" "str q1, [x9, #80]\n" "str q2, [x9, #96]\n" "str q3, [x9, #112]\n"
        "str q0, [x9, #128]\n" "str q1, [x9, #144]\n" "str q2, [x9, #160]\n" "str q3, [x9, #176]\n"
        "str q0, [x9, #192]\n" "str q1, [x9, #208]\n" "str q2, [x9, #224]\n" "str q3, [x9, #240]\n"
        "subs w10, w10, #1\n" "b.ne 1b\n"
        : [o] "+w"(acc) : [b] "r"(buf), [n] "r"(N)
        : "v0","v1","v2","v3","x9","w10","memory","cc");
    out[0] = acc;
}

// pass2-like mix: 16 ldr.q + 16 str.q + 12 fadd + 8 fmul + 4 fmla (all independent-ish)
static void mix(double *out)
{
    static double buf[16384] __attribute__((aligned(64)));
    for (int i = 0; i < 16384; i++) buf[i] = i * 1e-9;
    double acc = 0;
    __asm__ volatile(
        "mov x9, %[b]\n" "add x11, x9, #4096\n" "mov w10, %w[n]\n"
        "movi v28.2d, #0\n" "movi v29.2d, #0\n"
        "1:\n"
        "ldr q0, [x9]\n" "ldr q1, [x9, #16]\n" "ldr q2, [x9, #32]\n" "ldr q3, [x9, #48]\n"
        "ldr q4, [x9, #64]\n" "ldr q5, [x9, #80]\n" "ldr q6, [x9, #96]\n" "ldr q7, [x9, #112]\n"
        "ldr q8, [x9, #128]\n" "ldr q9, [x9, #144]\n" "ldr q10, [x9, #160]\n" "ldr q11, [x9, #176]\n"
        "ldr q12, [x9, #192]\n" "ldr q13, [x9, #208]\n" "ldr q14, [x9, #224]\n" "ldr q15, [x9, #240]\n"
        "fadd v16.2d, v0.2d, v1.2d\n" "fsub v17.2d, v2.2d, v3.2d\n"
        "fadd v18.2d, v4.2d, v5.2d\n" "fsub v19.2d, v6.2d, v7.2d\n"
        "fadd v20.2d, v8.2d, v9.2d\n" "fsub v21.2d, v10.2d, v11.2d\n"
        "fadd v22.2d, v12.2d, v13.2d\n" "fsub v23.2d, v14.2d, v15.2d\n"
        "fmul v24.2d, v16.2d, v17.2d\n" "fmla v24.2d, v18.2d, v19.2d\n"
        "fmul v25.2d, v20.2d, v21.2d\n" "fmla v25.2d, v22.2d, v23.2d\n"
        "str q16, [x11]\n" "str q17, [x11, #16]\n" "str q18, [x11, #32]\n" "str q19, [x11, #48]\n"
        "str q20, [x11, #64]\n" "str q21, [x11, #80]\n" "str q22, [x11, #96]\n" "str q23, [x11, #112]\n"
        "str q24, [x11, #128]\n" "str q25, [x11, #144]\n" "str q26, [x11, #160]\n" "str q27, [x11, #176]\n"
        "str q16, [x11, #192]\n" "str q17, [x11, #208]\n" "str q18, [x11, #224]\n" "str q19, [x11, #240]\n"
        "fadd v28.2d, v28.2d, v24.2d\n" "fadd v29.2d, v29.2d, v25.2d\n"
        "subs w10, w10, #1\n" "b.ne 1b\n"
        "fadd v28.2d, v28.2d, v29.2d\n" "fmov %d[o], d28\n"
        : [o] "+w"(acc) : [b] "r"(buf), [n] "r"(N/16)
        : "v0","v1","v2","v3","v4","v5","v6","v7","v8","v9","v10","v11",
          "v12","v13","v14","v15","v16","v17","v18","v19","v20","v21","v22","v23",
          "v24","v25","v28","v29","x9","x11","w10","cc");
    out[0] = acc;
}

// streaming read over working set with 4 independent accumulators
static void t_stream(double *buf, size_t nbytes, int reps, double *out)
{
    size_t n = nbytes / 16;
    double acc = 0;
    for (int r = 0; r < reps; r++) {
        __asm__ volatile(
            "mov x9, %[b]\n" "mov x10, %[n]\n"
            "movi v28.2d, #0\n" "movi v29.2d, #0\n" "movi v30.2d, #0\n" "movi v31.2d, #0\n"
            "1:\n"
            "ldp q0, q1, [x9], #32\n"
            "ldp q2, q3, [x9], #32\n"
            "ldp q4, q5, [x9], #32\n"
            "ldp q6, q7, [x9], #32\n"
            "fadd v28.2d, v28.2d, v0.2d\n" "fadd v29.2d, v29.2d, v2.2d\n"
            "fadd v30.2d, v30.2d, v4.2d\n" "fadd v31.2d, v31.2d, v6.2d\n"
            "fadd v28.2d, v28.2d, v1.2d\n" "fadd v29.2d, v29.2d, v3.2d\n"
            "fadd v30.2d, v30.2d, v5.2d\n" "fadd v31.2d, v31.2d, v7.2d\n"
            "subs x10, x10, #8\n" "b.ne 1b\n"
            "fadd v28.2d, v28.2d, v29.2d\n" "fadd v30.2d, v30.2d, v31.2d\n"
            "fadd v28.2d, v28.2d, v30.2d\n" "fmov %d[o], d28\n"
            : [o] "+w"(acc) : [b] "r"(buf), [n] "r"(n)
            : "v0","v1","v2","v3","v4","v5","v6","v7","v28","v29","v30","v31","x9","x10","cc");
    }
    out[0] = acc;
}

static void t_copy(double *dst, double *src, size_t nbytes, int reps, double *out)
{
    size_t n = nbytes / 16;
    for (int r = 0; r < reps; r++) {
        __asm__ volatile(
            "mov x9, %[s]\n" "mov x11, %[d]\n" "mov x10, %[n]\n"
            "1:\n"
            "ldp q0, q1, [x9], #32\n" "ldp q2, q3, [x9], #32\n"
            "stp q0, q1, [x11], #32\n" "stp q2, q3, [x11], #32\n"
            "subs x10, x10, #4\n" "b.ne 1b\n"
            : : [s] "r"(src), [d] "r"(dst), [n] "r"(n)
            : "v0","v1","v2","v3","x9","x10","x11","memory","cc");
    }
    out[0] = 0;
}

int main(void)
{
    double out[2], t;

    printf("%-10s %9s %9s %9s\n", "test", "ms", "ns/op", "cyc/op");
#define RUN(name, label, ops)                                            \
    t = now_sec(); name(out); t = now_sec() - t;                          \
    printf("%-10s %9.3f %9.4f %9.3f\n", label, t*1e3, t*1e9/(ops),      \
           t*GHZ*1e9/(ops));

    RUN(lat_fadd, "lat_fadd", (double)N) RUN(lat_fmul, "lat_fmul", (double)N)
    RUN(lat_fmla, "lat_fmla", (double)N)
    RUN(tp_fadd, "tp_fadd", (double)N*8) RUN(tp_fmul, "tp_fmul", (double)N*8)
    RUN(tp_fmla, "tp_fmla", (double)N*8)
    RUN(ld_q, "ld_q", (double)N*16)
    RUN(st_q, "st_q", (double)N*16)
    RUN(mix, "mix", (double)N)
    RUN(tp_fadd, "tp_fadd_b", (double)N*8)

    struct { const char *name; size_t bytes; int reps; } lv[] = {
        {"L1-8K", 8*1024, 800}, {"L1-32K", 32*1024, 250}, {"L2-128K", 128*1024, 64},
        {"L2-256K", 256*1024, 32}, {"L2-512K", 512*1024, 16}, {"DRAM-16M", 16*1024*1024, 1},
    };
    double *b = aligned_alloc(64, 16*1024*1024 + 64);
    double *d = aligned_alloc(64, 16*1024*1024 + 64);
    for (size_t i = 0; i < 16*1024*1024/8; i++) b[i] = i * 1e-9;
    for (int i = 0; i < (int)(sizeof lv/sizeof lv[0]); i++) {
        double r;
        t = now_sec(); t_stream(b, lv[i].bytes, lv[i].reps, &r); t = now_sec() - t;
        double tot = (double)lv[i].bytes * lv[i].reps;
        printf("read %-9s %9.3f ms %6.2f GB/s  %5.2f B/cyc\n", lv[i].name, t*1e3, tot/t/1e9, tot/t/GHZ/1e9);
        t = now_sec(); t_copy(d, b, lv[i].bytes, lv[i].reps, &r); t = now_sec() - t;
        printf("copy %-9s %9.3f ms %6.2f GB/s  %5.2f B/cyc(rw)\n", lv[i].name, t*1e3, 2*tot/t/1e9, 2*tot/t/GHZ/1e9);
    }
    free(b); free(d);
    return 0;
}
