// ubench.c -- Cortex-A53 microbenchmark: DP FP pipeline throughput + memory levels.
// Standalone static aarch64 binary (built with the same zig cc -mcpu=cortex_a53).
// Usage: ubench
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

#define N 2000000

// ---- 1. pure FP throughput: 8 independent 2-double chains ----
#define FPTEST(name, insn, mac)                                          \
    static void name(double *out)                                        \
    {                                                                    \
        double r0[2] = {1.0, 2.0}, r1[2] = {1.1, 2.1}, r2[2] = {1.2, 2.2}, r3[2] = {1.3, 2.3}; \
        double r4[2] = {1.4, 2.4}, r5[2] = {1.5, 2.5}, r6[2] = {1.6, 2.6}, r7[2] = {1.7, 2.7}; \
        double c[2] = {1e-9, 1e-9};                                      \
        for (int i = 0; i < N; i++) {                                    \
            __asm__ volatile(                                            \
                "0:\n"                                                   \
                insn " v0.2d, v0.2d, v16.2d\n"                           \
                insn " v1.2d, v1.2d, v16.2d\n"                           \
                insn " v2.2d, v2.2d, v16.2d\n"                           \
                insn " v3.2d, v3.2d, v16.2d\n"                           \
                insn " v4.2d, v4.2d, v16.2d\n"                           \
                insn " v5.2d, v5.2d, v16.2d\n"                           \
                insn " v6.2d, v6.2d, v16.2d\n"                           \
                insn " v7.2d, v7.2d, v16.2d\n"                           \
                "subs %w8, %w8, #1\n"                                    \
                "b.ne 0b\n"                                              \
                : [cnt] "+r"(i)                                          \
                : "w"(*r0), "w"(*r1), "w"(*r2), "w"(*r3),                \
                  "w"(*r4), "w"(*r5), "w"(*r6), "w"(*r7), "w"(*c)        \
                : "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "cc"); \
        }                                                                \
        out[0] = r0[0] + r1[0] + r2[0] + r3[0];                          \
        out[1] = r4[0] + r5[0] + r6[0] + r7[0];                          \
    }

// The asm clobbers r0..r7 but the C variables must be fed in; use a simpler
// approach: pass the constant via v16 and let objects stay in registers.
#define FPTEST2(name, insn)                                              \
    static void name(double *out)                                        \
    {                                                                    \
        double acc = 1.0, c = 1e-9;                                      \
        __asm__ volatile(                                                \
            "fmov d16, %d[cin]\n"                                        \
            "fmov d0, %d[acc]\n"                                         \
            "fmov d1, d0\n" "fmov d2, d0\n" "fmov d3, d0\n"              \
            "fmov d4, d0\n" "fmov d5, d0\n" "fmov d6, d0\n" "fmov d7, d0\n" \
            "mov w9, %w[n]\n"                                            \
            "1:\n"                                                       \
            insn " v0.2d, v0.2d, v16.2d\n"                               \
            insn " v1.2d, v1.2d, v16.2d\n"                               \
            insn " v2.2d, v2.2d, v16.2d\n"                               \
            insn " v3.2d, v3.2d, v16.2d\n"                               \
            insn " v4.2d, v4.2d, v16.2d\n"                               \
            insn " v5.2d, v5.2d, v16.2d\n"                               \
            insn " v6.2d, v6.2d, v16.2d\n"                               \
            insn " v7.2d, v7.2d, v16.2d\n"                               \
            "subs w9, w9, #1\n"                                          \
            "b.ne 1b\n"                                                  \
            "fadd d0, d0, d1\n"                                          \
            "fadd d2, d2, d3\n"                                          \
            "fadd d4, d4, d5\n"                                          \
            "fadd d6, d6, d7\n"                                          \
            "fadd d0, d0, d2\n"                                          \
            "fadd d4, d4, d6\n"                                          \
            "fadd d0, d0, d4\n"                                          \
            "fmov %d[o], d0\n"                                           \
            : [o] "+w"(acc)                                              \
            : [acc] "w"(acc), [cin] "w"(c), [n] "r"(N)                   \
            : "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v16", "w9", "cc"); \
        out[0] = acc;                                                    \
    }

FPTEST2(t_fadd, "fadd")
FPTEST2(t_fmul, "fmul")
FPTEST2(t_fmla, "fmla")

// ---- 2. load throughput (L1 hit): 16 independent LDR Q per iter ----
static void t_ldrq(double *out)
{
    static double buf[4096] __attribute__((aligned(64)));
    for (int i = 0; i < 4096; i++) buf[i] = i * 1e-9;
    double acc = 0;
    __asm__ volatile(
        "mov x9, %[b]\n"
        "mov w10, %w[n]\n"
        "movi v31.2d, #0\n"
        "1:\n"
        "ldr q0, [x9]\n" "ldr q1, [x9, #16]\n" "ldr q2, [x9, #32]\n" "ldr q3, [x9, #48]\n"
        "ldr q4, [x9, #64]\n" "ldr q5, [x9, #80]\n" "ldr q6, [x9, #96]\n" "ldr q7, [x9, #112]\n"
        "ldr q8, [x9, #128]\n" "ldr q9, [x9, #144]\n" "ldr q10, [x9, #160]\n" "ldr q11, [x9, #176]\n"
        "ldr q12, [x9, #192]\n" "ldr q13, [x9, #208]\n" "ldr q14, [x9, #224]\n" "ldr q15, [x9, #240]\n"
        "fadd v31.2d, v31.2d, v0.2d\n" "fadd v31.2d, v31.2d, v1.2d\n"
        "fadd v31.2d, v31.2d, v2.2d\n" "fadd v31.2d, v31.2d, v3.2d\n"
        "fadd v31.2d, v31.2d, v4.2d\n" "fadd v31.2d, v31.2d, v5.2d\n"
        "fadd v31.2d, v31.2d, v6.2d\n" "fadd v31.2d, v31.2d, v7.2d\n"
        "fadd v31.2d, v31.2d, v8.2d\n" "fadd v31.2d, v31.2d, v9.2d\n"
        "fadd v31.2d, v31.2d, v10.2d\n" "fadd v31.2d, v31.2d, v11.2d\n"
        "fadd v31.2d, v31.2d, v12.2d\n" "fadd v31.2d, v31.2d, v13.2d\n"
        "fadd v31.2d, v31.2d, v14.2d\n" "fadd v31.2d, v31.2d, v15.2d\n"
        "subs w10, w10, #1\n"
        "b.ne 1b\n"
        "fmov %d[o], d31\n"
        : [o] "+w"(acc)
        : [b] "r"(buf), [n] "r"(N)
        : "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10", "v11",
          "v12", "v13", "v14", "v15", "v31", "x9", "w10", "cc");
    out[0] = acc;
}

// ---- 3. LDP Q throughput ----
static void t_ldpq(double *out)
{
    static double buf[4096] __attribute__((aligned(64)));
    for (int i = 0; i < 4096; i++) buf[i] = i * 1e-9;
    double acc = 0;
    __asm__ volatile(
        "mov x9, %[b]\n"
        "mov w10, %w[n]\n"
        "movi v31.2d, #0\n"
        "1:\n"
        "ldp q0, q1, [x9]\n" "ldp q2, q3, [x9, #32]\n"
        "ldp q4, q5, [x9, #64]\n" "ldp q6, q7, [x9, #96]\n"
        "ldp q8, q9, [x9, #128]\n" "ldp q10, q11, [x9, #160]\n"
        "ldp q12, q13, [x9, #192]\n" "ldp q14, q15, [x9, #224]\n"
        "fadd v31.2d, v31.2d, v0.2d\n" "fadd v31.2d, v31.2d, v2.2d\n"
        "fadd v31.2d, v31.2d, v4.2d\n" "fadd v31.2d, v31.2d, v6.2d\n"
        "fadd v31.2d, v31.2d, v8.2d\n" "fadd v31.2d, v31.2d, v10.2d\n"
        "fadd v31.2d, v31.2d, v12.2d\n" "fadd v31.2d, v31.2d, v14.2d\n"
        "subs w10, w10, #1\n"
        "b.ne 1b\n"
        "fmov %d[o], d31\n"
        : [o] "+w"(acc)
        : [b] "r"(buf), [n] "r"(N)
        : "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10", "v11",
          "v12", "v13", "v14", "v15", "v31", "x9", "w10", "cc");
    out[0] = acc;
}

// ---- 4. FP+load mix: mimic pass2 ratio (62 fp : 16 ld : 16 st) ----
// 16 ldr + 16 str + 8 fadd + 8 fsub + 8 fmul + 8 fmla per iter
static void t_mix(double *out)
{
    static double buf[8192] __attribute__((aligned(64)));
    for (int i = 0; i < 8192; i++) buf[i] = i * 1e-9;
    double acc = 0;
    __asm__ volatile(
        "mov x9, %[b]\n"
        "add x11, x9, #2048\n"
        "mov w10, %w[n]\n"
        "movi v31.2d, #0\n"
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
        "fadd v26.2d, v24.2d, v25.2d\n" "fsub v27.2d, v24.2d, v25.2d\n"
        "fmul v28.2d, v26.2d, v27.2d\n" "fmla v28.2d, v16.2d, v20.2d\n"
        "str q16, [x11]\n" "str q17, [x11, #16]\n" "str q18, [x11, #32]\n" "str q19, [x11, #48]\n"
        "str q20, [x11, #64]\n" "str q21, [x11, #80]\n" "str q22, [x11, #96]\n" "str q23, [x11, #112]\n"
        "str q24, [x11, #128]\n" "str q25, [x11, #144]\n" "str q26, [x11, #160]\n" "str q27, [x11, #176]\n"
        "str q28, [x11, #192]\n" "str q16, [x11, #208]\n" "str q17, [x11, #224]\n" "str q18, [x11, #240]\n"
        "fadd v31.2d, v31.2d, v28.2d\n"
        "subs w10, w10, #1\n"
        "b.ne 1b\n"
        "fmov %d[o], d31\n"
        : [o] "+w"(acc)
        : [b] "r"(buf), [n] "r"(N / 16)
        : "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v8", "v9", "v10", "v11",
          "v12", "v13", "v14", "v15", "v16", "v17", "v18", "v19", "v20", "v21", "v22",
          "v23", "v24", "v25", "v26", "v27", "v28", "v31", "x9", "x11", "w10", "cc");
    out[0] = acc;
}

// ---- 5. streaming read bandwidth over a working set of size S (sum of 2-double adds) ----
static void t_stream(double *buf, size_t nbytes, int reps, double *out)
{
    size_t n = nbytes / 16;   // Q words
    double acc = 0;
    for (int r = 0; r < reps; r++) {
        __asm__ volatile(
            "mov x9, %[b]\n"
            "mov x10, %[n]\n"
            "movi v31.2d, #0\n"
            "1:\n"
            "ldp q0, q1, [x9], #32\n"
            "ldp q2, q3, [x9], #32\n"
            "ldp q4, q5, [x9], #32\n"
            "ldp q6, q7, [x9], #32\n"
            "fadd v31.2d, v31.2d, v0.2d\n"
            "fadd v31.2d, v31.2d, v2.2d\n"
            "fadd v31.2d, v31.2d, v4.2d\n"
            "fadd v31.2d, v31.2d, v6.2d\n"
            "subs x10, x10, #8\n"
            "b.ne 1b\n"
            "fmov %d[o], d31\n"
            : [o] "+w"(acc)
            : [b] "r"(buf), [n] "r"(n)
            : "v0", "v1", "v2", "v3", "v4", "v5", "v6", "v7", "v31", "x9", "x10", "cc");
    }
    out[0] = acc;
}

// ---- 6. copy bandwidth (ldp + stp) over working set ----
static void t_copy(double *dst, double *src, size_t nbytes, int reps, double *out)
{
    size_t n = nbytes / 16;
    for (int r = 0; r < reps; r++) {
        __asm__ volatile(
            "mov x9, %[s]\n"
            "mov x11, %[d]\n"
            "mov x10, %[n]\n"
            "1:\n"
            "ldp q0, q1, [x9], #32\n"
            "ldp q2, q3, [x9], #32\n"
            "stp q0, q1, [x11], #32\n"
            "stp q2, q3, [x11], #32\n"
            "subs x10, x10, #4\n"
            "b.ne 1b\n"
            :
            : [s] "r"(src), [d] "r"(dst), [n] "r"(n)
            : "v0", "v1", "v2", "v3", "x9", "x10", "x11", "memory", "cc");
    }
    out[0] = 0;
}

int main(void)
{
    double out[2];
    double t;
    printf("clock=%lu Hz reported\n", (unsigned long)CLOCKS_PER_SEC);

#define RUN(name, label, ops)                                            \
    t = now_sec(); name(out); t = now_sec() - t;                          \
    printf("%-12s %8.3f ms  %6.3f ns/op  %6.2f ops/cycle@2016MHz\n",      \
           label, t * 1e3, t * 1e9 / (ops), (ops) / (t * 2.016e9));

    RUN(t_fadd, "fadd.2d", (double)N * 8)
    RUN(t_fmul, "fmul.2d", (double)N * 8)
    RUN(t_fmla, "fmla.2d", (double)N * 8)
    RUN(t_ldrq, "ldr.q", (double)N * 16)
    RUN(t_ldpq, "ldp.q", (double)N * 8)
    RUN(t_mix, "mix", (double)(N / 16))

    // memory levels
    struct { const char *name; size_t bytes; int reps; } lv[] = {
        {"L1-16K", 16 * 1024, 400},
        {"L1-32K", 32 * 1024, 200},
        {"L2-128K", 128 * 1024, 60},
        {"L2-256K", 256 * 1024, 30},
        {"L2-512K", 512 * 1024, 16},
        {"DRAM-16M", 16 * 1024 * 1024, 1},
    };
    double *b = aligned_alloc(64, 16 * 1024 * 1024 + 64);
    double *d = aligned_alloc(64, 16 * 1024 * 1024 + 64);
    for (size_t i = 0; i < 16 * 1024 * 1024 / 8; i++) b[i] = i * 1e-9;
    for (int i = 0; i < (int)(sizeof lv / sizeof lv[0]); i++) {
        double r;
        t = now_sec();
        t_stream(b, lv[i].bytes, lv[i].reps, &r);
        t = now_sec() - t;
        double tot = (double)lv[i].bytes * lv[i].reps;
        printf("read %-10s %8.3f ms  %6.2f GB/s  %5.2f B/cycle\n",
               lv[i].name, t * 1e3, tot / t / 1e9, tot / t / 2.016e9);
        t = now_sec();
        t_copy(d, b, lv[i].bytes, lv[i].reps, &r);
        t = now_sec() - t;
        printf("copy %-10s %8.3f ms  %6.2f GB/s(+%s) %5.2f B/cycle rw\n",
               lv[i].name, t * 1e3, 2 * tot / t / 1e9, lv[i].name, 2 * tot / t / 2.016e9);
    }
    free(b); free(d);
    return 0;
}
