# ARM64 vs x86 performance of the asm→C gwnum port

Goal: separate *port/code quality* from *hardware differences*. All numbers were
measured with one small benchmark (`_work/bench/gwnum-bench.c`, public gwnum API:
`gwmul3` of two random numbers, single thread) linked against:

* **our asm→C converted gwnum** — the library used by the released ARM64 binary
  (aarch64: `-O3 -ffp-contract=on -mcpu=cortex_a53`; x86-64: MinGW `-O2`), and
* the **official x86-64 asm gwnum** shipped with the Prime95 source
  (`gwnum/gwnum64.lib`, MSVC 2022 build, linked into the same benchmark).

On x86 the CPU profile can be forced to the 128-bit-only ARM64 profile
(`CPU_FLAGS=0x3F2F`, `CPU_ARCHITECTURE_OTHER`) so that both builds select the
same FFT layout (`type-1 FFT length 16K, Pass1=64, Pass2=256, clm=4`). The
official lib was measured both with that forced 128-bit profile (apples-to-
apples with the converted C code) and with its native ISA (FMA3/AVX2 on this
machine).

## Raw results

p = 264,169 (16K FFT) and p = 7,831,391 (400K FFT), one multiply per iteration,
single thread. Cycles: `QueryThreadCycleTime` on x86 (real core cycles, turbo
included), `time × 2016 MHz` on the RK3528 (fixed-frequency A53), `time ×
nominal clock` for the phone (its `scaling_cur_freq` sampling was unreliable
during the short runs, so those two rows are marked ~).

| platform | core / clock | 16K µs/op | 16K cycles/op | 400K µs/op | 400K cycles/op |
|---|---|---|---|---|---|
| x86, official asm, native ISA (FMA3/AVX2) | Arrow Lake, ~3.6 GHz TSC (turbo) | 23.9 | 86,092 | 845 | 3,045,263 |
| x86, official asm, forced 128-bit (SSE2) | same | 43.4 | 158,052 | 1,533 | 5,555,738 |
| **x86, our converted C, forced 128-bit** | same | **42.4–45.2** | **154,200–156,000** | **1,515** | **5,501,705** |
| ARM, our converted C | Cortex-X925, ~3.5 GHz | 55.7 | ~195,000 | — | — |
| ARM, our converted C | Cortex-A725, 2.1 GHz | 160.0 | ~336,000 | — | — |
| ARM, our converted C | Cortex-A53, 2.016 GHz | 1,147 | 2,312,468 | 36,463 | 73,509,562 |

## What the numbers say

1. **The asm→C conversion costs ~nothing.** On the same x86 machine, with the
   same 128-bit ISA profile, the converted C kernels are *as fast as the original
   hand-written asm*: 154–156k vs 158k cycles at 16K (C is ~1% faster) and
   5.50M vs 5.56M at 400K. So the port's code is not the bottleneck — the ARM
   version's code is essentially "the x86 version" translated to NEON, at
   parity per operation.

2. **The x86 version's real speed edge is ISA width + clock.** Its native build
   uses 256-bit FMA3 vectors instead of 128-bit SSE2/NEON: 23.9 vs 43.4 µs at
   16K (1.8×), 845 vs 1,533 µs at 400K (1.8×). Our ARM port uses 128-bit NEON
   (the A53 and the farm box have no SVE); a SVE/SVE2 path (Cortex-X925 class)
   is the obvious future optimization.

3. **Per-clock, per-core** (dividing the cycle counts):
   * Cortex-X925 ≈ **80%** of the Arrow Lake core — modern ARM big cores are
     close to x86 cores at the same clock;
   * Cortex-A725 ≈ 47%;
   * Cortex-A53 (the RK3528 farm box) ≈ **7%** — an in-order, 2-wide core,
     i.e. a hardware class difference, not a port issue.

4. **With real clocks** (A53 2.016 GHz vs x86 ~4.5–5.5 GHz boost): the box's
   single core is roughly **1/25** of one x86 core on this workload; a phone
   X925 (~3.5 GHz) is about 0.6× one x86 core. The farm's advantage is
   cheap 24/7 operation, not per-core speed.

Reproduce: `_work/bench/gwnum-bench.c` (x86 MSVC: `/O2` + `gwnum64.lib` +
`polymult64.lib` + `advapi32.lib`; x86 MinGW: converted-C libs; aarch64:
`zig cc -target aarch64-linux-musl` + `lib-a53/gwnum.a`).
