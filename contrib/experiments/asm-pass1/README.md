# contrib/experiments/asm-pass1 — hand-written NEON for the pass-1 copy loops

**Verdict: technically correct, measurably useless end-to-end (≈0%, within noise).
Do not repeat this shape of work; see "What this tells us" below.**

This directory archives a completed experiment (2026-11-16) on the RK3528
ARM64 farm, so that the measurement tools survive and the negative result is
not rediscovered the hard way.

## What was done

`r4pass1_asm.S` replaces the two pass-1 scratch↔FFT cache-line swizzle loops
`copy_scratch_data_to_fft` / `copy_fft_data_to_scratch` (macros in `xmult.h`,
instantiated in ~30 pass-1 variants) for the `clm=4` kernels (16K/48K/80K np)
with a software-pipelined AArch64 NEON kernel (16 loads issued before 16
stores, v0–v15, post-indexed addressing). Gated by
`-DGWNUM_ASM_PASS1_COPY`; the C version stays the default.

## Measured results (RK3528, 4x Cortex-A53 @ 2016 MHz, farm paused)

| measurement | baseline C | asm | ratio |
|---|---|---|---|
| isolated copy kernel (`copytest` microbench, same geometry) | 1541 ns/call | 1131 ns/call | **1.36x** |
| end-to-end multiply, 16K FFT (exp 264169), interleaved A/B | 1.0775–1.0825 ms | 1.0898–1.0951 ms | **≈0 (noise)** |
| end-to-end multiply, 400K FFT (exp 7831391) | 35.69–36.14 ms | untouched | — |

Bit-exactness: the asm build produced identical digests to the baseline on both
sizes (40+ runs), and the differential unit test `copytest` passed for two
geometries with 1 MB random buffers.

## Why the 1.36x did not survive

- The ported loops are only ~14% of pass-1 samples and ~3–7% of a full
  multiply, so the theoretical ceiling was ≈1%.
- They are **bandwidth-bound, not scheduling-bound**: 8192 bytes in ~1.13 µs =
  7.4 GB/s effective against an 8320-byte-stride destination. A hand-written
  copy cannot do better than the 1.36x measured, and that gain is smaller than
  the run-to-run noise plus code-layout effects (the binary shrank by 158 KB
  because the copy loops are instantiated in ~30 variants, so the layout of all
  pass-1 hot loops moved).

## Corrected hot-spot profile (this experiment's most useful output)

An earlier top-50-PC list suggested 16K = 72.5% pass1 / 25.6% xi2 / 1.9%
pass2. That was an artifact: pass-2's samples are spread over hundreds of
distinct PCs in its 4 KB unrolled body, so none of them reached a top-50 list.
Attributing **every** sample (1305 samples, SIGPROF 1 ms):

| 16K (M264169) | % |  | 400K (M7831391) | % |
|---|---|---|---|---|
| `xfft_r4_16K_np_8_4_CORE_i` | **47.4** |  | `xfft_r4dwpn_400K_1280_4_CORE_i` | **51.6** |
| `r4_pass2_8_levels_complex` | **33.3** |  | `r4_pass2_1280_levels_complex` | **38.7** |
| `xi2CORE` | 14.7 |  | `xi3CORE` | 8.8 |
| `xpass2_r4_8_levels` | 3.1 |  | `xpass2_r4_1280_levels` | 0.7 |

Loop map of the 16K pass-1 function (5428 B, 1357 instrs, 618 samples): the
hottest single block is only 11% of total runtime, the code is flat and
memory/stall bound (IPC ≈ 0.3–0.5). That is why porting any one sub-macro
cannot pay off.

## What this tells us

1. **Pass-2 is the largest untapped target** (33% at 16K, 39% at 400K) — a
   pass-2 rewrite is the right next experiment (in progress separately).
2. Sub-macro-level asm on this codebase has a **1–2% ceiling**, below device
   noise. Don't do it again without a plan that covers a much larger fraction.
3. The real limiter is the L1 working set: 128 KB of FFT data + 128 KB of
   scratch vs 32 KB L1D. Beating it needs an **algorithm/layout change**
   (blocking), not instruction scheduling.
4. Method discipline that must be kept: farm paused with `kill -STOP` /
   `kill -CONT` (always), `taskset 4` bitmask, interleaved A/B, medians, and
   bit-exact digests with a **fixed-length** hash (hashing `|sign|` limbs of a
   freshly-allocated giant hashes uninitialised memory and is
   non-deterministic — that mistake cost a day here).

## Archived tools (reusable)

- `fftbench_asm.c` + `prof_asm.c` — deterministic digest benchmark / PC-sampling
  profiler. `fftbench_asm` prints `timeA_ms` (x=3.0 chain), `timeB_ms`
  (full-size pseudo-random operand) and a digest over the raw FFT doubles and
  giant limbs; its digest is stable across runs and identical between builds
  that compute the same values.
- `np_*.py` — decode sampled PCs into per-symbol / per-loop / per-region
  profiles (no binutils for aarch64 needed; ELF symtab is parsed directly).
- `copytest.c` — differential test for a data-movement replacement.
- `build_base.cmd` / `build_asm.cmd` / `build_bench.cmd` — reference build
  recipes (zig cc cross-compile, forward-slash OBJDIR/LIBDIR).
- `ab_bench.sh` / `ab2.sh` / `run_asm.sh` / `prof_run.sh` — on-device
  measurement scripts (farm pause/resume protocol included).
- `r4pass1_asm.S` — the asm kernel itself (kept for reference; not wired into
  the default build).

Correctness gate for any such experiment: `../validation/` (three layers).
