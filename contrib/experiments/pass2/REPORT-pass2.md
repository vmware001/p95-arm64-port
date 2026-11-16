# Pass-2 FFT kernel optimisation experiment (16K / 400K, ARM64 Cortex-A53)

Owner: pass2-opt agent. Workspace: `D:\data\code\p95\_work\pass2\` only.
Git tree `p95v3107b02.arm64.source` untouched (`git status --porcelain` empty).
Production binary `_work/mprime/mprime-mc4`, `_work/bench/`, `_work/gwnum/lib-*`,
`_work/asm/` untouched. Farm left running (mprime alive, 4 workers).

## Verdict

**No measurable speedup found.** Three legal, bit-exact data-movement changes were
built and measured end-to-end, interleaved, on the box:

* 4-iteration inner butterfly loops unrolled x2 (more ILP): **+0.5% slower**
* restored twiddle software prefetch (`prfm pldl1keep`, 256 B lead) in place of
  the mechanically-converted no-op `xxprefetch`: **+0.4% slower**
* `__restrict` on the four butterfly pointers (to let the compiler hoist the 8
  interleaved loads above the 8 stores): compiled to **byte-identical code**
  (same MD5 for the whole 4096-byte function) - clang already proved the four
  64-byte regions disjoint, so the interleaved schedule is its deliberate choice,
  not an aliasing artifact. Not worth a timing run.

All are within/near the noise floor but consistently *not* faster. The
characterisation below explains why: the kernel is already running at the rate
this in-order core delivers for its instruction mix, and there is no ILP budget
left (32 NEON registers, each butterfly needs ~30) to hide anything.

## Characterisation

PC-sampling profile (1 ms SIGPROF; 5961 samples over 18000 full-size 16K
multiplies, exp 264169):

| function | share |
|---|---|
| `xfft_r4_16K_np_8_4_CORE_i` (pass1) | 49.7% |
| **`r4_pass2_8_levels_complex`** (pass2) | **33.7%** |
| `xi2CORE` | 12.7% |
| `xpass2_r4_8_levels` | 2.9% |
| `xgw_carries` | 0.3% |

(agrees with the previous agent's 47.4 / 33.3 / 14.7 / 3.1 profile)

`r4_pass2_8_levels_complex`: 4096 bytes, 1024 instructions, 8 executed loops
(+3 dead `with_mult`/`with_mulf` loops that a squaring never reaches).

Hot inner butterfly loop (`r4_x4cl_four_complex_djbfft` / `djbunfft`):
103 instructions = **62 DP vector FP ops** (16 fadd, 16 fsub, 12 fmul, 12 fmla,
6 fneg) + **34 memory ops** (16 `ldr q`, 16 `str q`, 2 address) + 1 branch.
No spills, no scalar fallback; all 32 NEON registers live. Samples are spread
almost uniformly across the FP and memory instructions (largest single hotspot:
a loop-entry load pair at 12.6% of its loop) - i.e. throughput/stall limited
everywhere, not one bad dependency.

Per multiply: 16 blocks x 1024 doubles = 16384 doubles = the whole 16K FFT;
256 butterfly calls per block -> 4096 calls, ~179 cycles/call at the nominal
2016 MHz (IPC ~0.55). 254K DP vector ops in total.

Microbenchmarks on the same box (farm frozen, pinned to cpu2):

| measurement | value |
|---|---|
| DP `fadd`/`fmul`/`fmla` throughput | 0.634 ns/op each (one DP pipe; 1.28 cyc/op @2016 MHz, 1.0 cyc/op if real clock ~1.58 GHz) |
| DP latency | fadd 2.25 ns, fmul 2.30 ns, fmla 4.51 ns |
| 128-bit L1 store throughput | 1.28 cyc/store |
| L1 read / L2-128K read / DRAM read | 8.9 / 6.8 / 2.6 GB/s |
| DRAM copy | 1.29 GB/s |
| **replica of pass2 mix** (16x128-bit ld + 16 st + 14 DP, all L1 hits) | **2.08 cyc/instr** |

The replica's 2.08 cycles/instruction equals the real kernel's 2.06: the kernel
is already at the hardware's delivered rate for this instruction mix, with no
cache-miss contribution left to recover (working set per block is only 4 KB and
L1-resident after the first pass).

Ceiling: the FP-pipe limit for pass2 is ~325K cycles vs ~734K used (2.25x).
Perfect FP/memory overlap would cap pass2 at ~1.45x = ~11% end-to-end, but the
hardware does not deliver that overlap for this mix, and two butterflies cannot
be interleaved (needs ~60 vector registers; ISA has 32) - which is why the x2
unroll is neutral-to-worse. A real win would need fewer DP ops (radix-8 /
twiddle-count reduction), which changes rounding and is forbidden by the
bit-exactness gate.

## Measurements (median of interleaved runs, farm SIGSTOPped, `taskset 4` = cpu2)

16K FFT, exp 264169, 1000 gwmul3 iterations, 3 reps x 5 runs = 15 runs/spec:

| spec | timeA ms | timeB ms | digest |
|---|---|---|---|
| `bench_base` (production) | **1.0637** | 1.0675 | 9eb6e6aa156860ac |
| harness v0 (same source, in-binary control) | 1.0703 | 1.0708 | 9eb6e6aa156860ac |
| v1 = inner loops unrolled x2 | 1.0760 (+0.5% vs v0) | 1.0747 | 9eb6e6aa156860ac |
| v4 = twiddle `prfm` L1, lead 256 B | 1.0752 (+0.4% vs v0) | 1.0776 | 9eb6e6aa156860ac |

400K FFT, exp 7831391, 80 iterations, 2 reps x 3 runs = 6 runs/spec:

| spec | timeA ms | timeB ms | digest |
|---|---|---|---|
| `bench_base` (production) | **35.091** | 35.313 | ae78e849c76aef67 |
| harness v0 | 35.289 (+0.56%) | 35.388 | ae78e849c76aef67 |
| v1 = 1280-level twiddle `prfm` | 35.227 (+0.39%) | 35.325 | ae78e849c76aef67 |

Run-to-run spread within a spec: +-0.5-1.5%. The in-binary control (v0) carries
+0.6% over the standalone binary (dispatcher + code layout), so variants must be
compared against v0, not against `bench_base`.

## Bit-exactness

Every timed run's digest was checked (the binary prints it). 16K and 400K
digests are identical between the production binary and every valid variant
(60 + 18 runs). The v2/v3 variants built from the 'big'-loop unrolling are
**invalid experiments**, not evidence: my generator assumed a 1-line body for
the inverse-1,2 loop, which actually has 2 lines (`djbunfft(...)` +
`rdi += XMM_SCD2`), so the generated code dropped the twiddle advance and the
digest changed (477e71c520261530). They were discarded.

No optimisation was adopted, so the layer-1 (`exact_strict`) and layer-2
(`mprime_known_results.sh`) gates were not re-run: the production build is
byte-for-byte the previously validated one.

## Files

Created under `_work/pass2/` (gwnum tree is a private copy; source repo untouched):

* `gwnum/xbasics.h` - added `GWNUM_ARM_PREFETCH` block: turns the converted-to-
  no-op `xprefetch`/`xxprefetch`/`xprefetchw` back into real PRFM hints (off by
  default, so the normal build is unchanged).
* `gwnum/r4/r4vwrap_8.c`, `r4vwrap_1280.c` - compile a variant of the pass-2
  source with renamed entry points (`-DVARIANT=vN -DSRC8=...`).
* `gwnum/r4/r4vdisp.c` - `GWNUM_P2VAR` env dispatcher, so all variants live in
  one binary and A/B comparisons share one code layout.
* `gwnum/r4/r4pass2_8_u_{loopc,all,big}.c` - generated restructured copies.
* `gwnum/r4/r4.h` - private copy with `__restrict` on the butterfly data
  pointers (variant v5; produced byte-identical code, experiment closed).
* `gen_var8.py`, `check_var.py`, `parse_ab.py`, `p2_analyze.py`.
* `ab2.sh`, `frozen.sh` - farm-frozen benchmark drivers (a hard watchdog resumes
  the farm after 420 s no matter what happens to the driver), `fftbench_asm.c` /
  `prof_asm.c` reused unmodified from `_work/asm`.
* `ubench.c`, `ubench2.c` - A53 microbenchmarks; `prof2.c` - hrtimer profiler.
* `build_base.cmd`, `build_var.cmd`, `build_restrict.cmd`, `build_misc.cmd`,
  `bin/`, `obj/`, `lib/`, `logs/`.

Build recipe: identical to production (`zig cc -target aarch64-linux-musl
-O3 -ffp-contract=on -mcpu=cortex_a53`, `Makefile.arm64`, OBJDIR/LIBDIR in the
private tree, polymult.c twice plain/+`-DSSE2`), then link the harness with the
variant objects before `libgwnum.a` so they override the stock members.

## Ceiling / what would be needed

Measured end-to-end potential is **< 1-2%** for any data-movement change. The
only lever that could move pass2 materially is reducing the DP operation count
(radix-8 butterflies, fewer twiddle multiplies), which changes rounding and is
therefore forbidden by the bit-exactness gate. A hand-written NEON port of the
hottest block is unlikely to beat clang here: the generated loop is already
optimal-shaped (62 FP ops, no spills, loads scheduled ahead of uses, `ldp`/`stp`
where the addresses allow, ~2.06 cycles/instruction = exactly what the A53
delivers for this mix as measured by an independent replica).

Caveat on absolute cycle figures: the nominal clock is 2016 MHz but the DP
throughput measures as 0.78 vector ops/cycle there, which suggests the real
clock may be ~1.58 GHz; all *ratios* above are clock-independent.

## Farm safety note

The bench driver freezes the farm with `SIGSTOP` and always sends `SIGCONT`.
Two incidents occurred while iterating on the harness and are disclosed for
transparency: (1) a 3-rep batch outlived the tool call that launched it and the
ssh client was killed before the trap ran, leaving mprime stopped for ~10 minutes
before I resumed it (`kill -CONT`); (2) one interrupted batch left its runs
partially measured while the farm was resumed - those runs were discarded and the
final data was re-taken with a single detached batch. The driver now runs
detached with a hard 420 s watchdog that resumes the farm no matter what. Final
state: mprime alive (pid 9216, 7 threads = 4 ECM workers), running at normal
~20 s per stage-1 block, no benchmark/scratch processes left on the box.

## Other observations worth passing on

* Pass1 (`xfft_r4_16K_np_8_4_CORE_i`) is 49.7% of the 16K multiply and has the
  same profile shape (uniform stalls, ~2 cycles/instruction), so the same
  "already at machine speed" conclusion probably applies there.
* The mechanical asm->C conversion disabled *all* the x86 prefetch macros and
  the whole `loops_init_prefetch` / `pfing` / `data_prefetch` machinery. On this
  in-order core that turned out not to matter for pass2 (working set 4 KB), but
  it may matter for kernels with larger streams (e.g. `xmult` / carries /
  400K pass1), which are outside this task.
* `xi2CORE` is 12.7% at 16K - the next-largest single target after pass1/pass2.
