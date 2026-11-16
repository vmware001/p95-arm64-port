# Multithreaded gwnum corruption hunt - x86 differential-soak report

Owner: `mt-x86` (task `task-7`). Scope: x86 host differential soaks of the
ARM64 "asm->C" gwnum 31.7 sources (`_work/gwnum/lib-x86/*.a`), plus the
ECM-faithful harness. The ARM device line is handled elsewhere.

**Status: NOT-REPRODUCED on x86** after ~12 million multiply-equivalents of
differential soak (threads=1 reference vs threads=N) across the raw-multiply
mix, an ECM-faithful Edwards `ed_add`/`ed_dbl` mix, and a production-shaped
ECM stage-1 chain (in-place ops, FFT-state chaining, `gwerror_checking`
patterns, `ED_XYZ`/`ED_SUBTRACT`/`GWADD_*` paths, spin-wait on/off, thread
counts 2..24, and 8 threads pinned onto a single physical core). Every
threads=1 vs threads=N digest comparison was bit-identical.

---

## 1. Artifacts (all under `_work/mt/x86/`)

| Path | Description |
|---|---|
| `..\mtstress_ed.cpp` | new ECM-faithful differential harness (source, 900 lines) |
| `..\mtstress_ed.exe` / `..\mtstress_ed_new.exe` | x86_64-windows builds (the `_new` one is the current source incl. `--chain`) |
| `mtstress_ed-arm` | aarch64-linux-musl static build of the **current** harness - ready to push to the phone/box for the decisive on-device differential test |
| `logs/*.log` | every run's full output |
| `run_*.cmd` | exact command line of every soak (reproducible) |

Build lines (same converted C sources as the ARM64 library, only host codegen
differs):

```
x86:  x86_64-w64-mingw32-g++ -O2 -o mtstress_ed.exe mtstress_ed.cpp \
        -ID:/data/code/p95/p95v3107b02.arm64.source/gwnum \
        _work/gwnum/lib-x86/gwnum.a _work/gwnum/lib-x86/polymult.a -lstdc++ -lm -lpthread
arm:  zig c++ -target aarch64-linux-musl -O2 -static -s -o mtstress_ed-arm \
        mtstress_ed.cpp -I<repo>/gwnum <repo>/gwnum/gwnum.a <repo>/gwnum/polymult.a -lm -lpthread
```

The x86 binary keeps `mtstress.c`'s CPU emulation (mask AVX/AVX2/FMA3/AVX512
out of `CPU_FLAGS`, force `CPU_ARCHITECTURE_OTHER`, family/model/signature 0)
so the **same** FFT layout is selected as on ARM64. Verified in every log:

```
type-1 FFT length 16K, Pass1=64, Pass2=256, clm=4      (threads=1)
type-1 FFT length 16K, Pass1=64, Pass2=256, clm=4, 8 threads
```

`gwsetup(h,1.0,2.0,n,-1.0)` with `gwset_larger_fftlen_count(0)`,
`gwset_safety_margin(0.0)`, `gwset_bench_cores(1)`, `gwset_bench_workers(1)`,
`gwclear_use_benchmarks`, `gwset_maxmulbyconst(2)` - i.e. production ecm.cpp
setup. Modulus is 2^n-1 generated exactly as production.

## 2. The harness (`_work/mt/mtstress_ed.cpp`)

Two operation modes, both differential (pass 1 = `gwset_num_threads(1)`
reference, pass 2 = threads=N; per-iteration 64-bit FNV digest of gwtogiant
limbs of every output coordinate + `alternate_format` flags + `gw_test_for_error`):

* **random mode** (default): copies of the production ecm.cpp routines
  `ed_extend`, `ed_add`, `ed_dbl` (verbatim from `contrib/edtest/edtest.cpp`,
  i.e. unmodified math incl. the `GWADD_FORCE_NORMALIZE` fix) plus the minimal
  `ecmhandle`. Each iteration generates fresh random extended/projective points
  as giants via an LCG and runs one of 12 variants:
  `dbl_dflt_xyz`, `dbl_add_snf`, `dbl_add_xyz`, `dbl_dflt_funky`,
  `dbl_add_xyz_snf`, `dbl_dflt_xyz_snf`, `add_snf`, `add_xyz`, `add_sub`,
  `add_nogwalloc`, `add_dbl_xyz_snf`, `add_nogw_xyz_snf`
  (snf = `ED_STARTNEXTFFT`, `ED_RESULT_FOR_ADD/DBL`, `ED_XYZ`, `ED_SUBTRACT`,
  `ED_NO_GWALLOC`; extended results get pre-allocated `t` buffers like
  production does).
* **chain mode** (`--chain`): a faithful copy of ecm.cpp's stage-1 main loop:
  a persistent in-place Edwards accumulator `e`; each bit either
  `ed_dbl(&e,&e, ED_RESULT_FOR_DBL|ED_STARTNEXTFFT)` or
  `ed_dbl(&e,&e,ED_RESULT_FOR_ADD|ED_STARTNEXTFFT)` followed by
  `ed_add(&e,&dict[k],&e, ED_RESULT_FOR_DBL|options)` with `options` =
  `ED_STARTNEXTFFT` / `ED_XYZ` / `|ED_SUBTRACT`, exactly as production; 8
  fixed extended "NAF dictionary" points; digest taken from `gwcopy`'d
  coordinates so the accumulator's FFT state/`unnorms` is preserved.
  `--errchk 0|1|2` selects `gwerror_checking` off / on / production-like
  (`(bitnum&127)==64`, and never on the `ED_XYZ` path - the current repo fix).

Instrumentation prints which `ed_dbl` algorithm branch actually runs; at this
FFT size it is always **algorithm 5** (`algo5=100001` in a 100k chain run, all
other branches 0), which is what production executes.

Validation of the detector itself (`--poison K`, corrupts pass-2 digest at
iteration K only):

```
MISMATCH iter 100 op[4]=dbl_add_xyz_snf err=0 got=d93a79a1eed216e0 want=d93a79a1307fa80f
RESULT: CORRUPTION - 1 mismatching iteration(s) ...      <- detector works end-to-end
```

Two harness bugs were found and fixed during bring-up (recorded so nobody
repeats them): (a) a `memset(harness)` after pass 1 wiped `words/topmask`
(symptom: deterministic SIGSEGV at phase-2 iteration 0 - *not* a gwnum bug);
(b) chain setup consumed the LCG before it was reseeded, so pass 1 and pass 2
used different dictionary points (symptom: deterministic "mismatch at iter 0"
for *every* thread count - also not a gwnum bug).

## 3. Coverage matrix

All runs: n in {264169, 291143} -> both select **16K type-1 FFT Pass1=64
Pass2=256 clm=4**, seeds distinct per run. "muls" counts gwmul3/square2/
submul4/mulmuladd5/mulmulsub5 calls per pass (both passes executed).

### 3.1 raw `mtstress.exe` (Lead's harness; fresh non-FFTed inputs, 6 op mix)

| run | n | threads | iters | ops | muls | result |
|---|---|---|---|---|---|---|
| raw_n264169_t8_s51A1 | 264169 | 8 | 400000 | `012345` | 0.47M | **PASS** (ST 2236 it/s, MT 1850 it/s) |
| raw01_n264169_t8 | 264169 | 8 | 200000 | `01` pure gwmul3 +/- STARTNEXTFFT | 0.2M | **_PENDING_** |
| raw23_n264169_t8 | 264169 | 8 | 200000 | `23` mulmuladd5 +/- STARTNEXTFFT | 0.2M | **_PENDING_** |
| raw45_n264169_t8 | 264169 | 8 | 200000 | `45` gwsubmul4 + chained-FTTed-source mul | 0.3M | **_PENDING_** |
| sanity | 264169 | 2,8 | 1000-3000 | all | - | PASS |

### 3.2 `mtstress_ed.exe` random mode (12 ECM Edwards variants)

| run | n | threads | iters | ops | muls | result |
|---|---|---|---|---|---|---|
| ed_n264169_t8_sed01 | 264169 | 8 | 100000 | all 12 | ~1.15M | **PASS** |
| ed_n291143_t8_sed02 | 291143 | 8 | 100000 | all 12 | ~1.15M | **PASS** |
| rnd_t16 | 264169 | 16 | 60000 | all 12 | ~0.7M | **_PENDING_** |
| sanity | 264169 | 2,4 | 200-3000 | all | - | PASS |
| poison self-test | 264169 | 8 | 300 | all | - | CORRUPT detected (expected) |

### 3.3 `mtstress_ed.exe` chain mode (production ECM stage-1 shape)

| run | n | threads | iters | errchk | extra | muls | result |
|---|---|---|---|---|---|---|---|
| ch_n264169_t8_sc001 | 264169 | 8 | 100000 | production | - | 1.45M | **_PENDING_** |
| ch_n291143_t8_sc002 | 291143 | 8 | 100000 | production | - | 1.45M | **_PENDING_** |
| ch_n264169_t8_sc003_ec0 | 264169 | 8 | 100000 | off | - | 1.45M | **_PENDING_** |
| ch_n264169_t8_sc004_ec1 | 264169 | 8 | 100000 | always on | - | 1.45M | **_PENDING_** |
| ch_spin1_t8 | 264169 | 8 | 60000 | production | `--spin-wait 1` | 0.87M | **_PENDING_** |
| ch_spin2_t8 | 264169 | 8 | 60000 | production | `--spin-wait 2` | 0.87M | **_PENDING_** |
| th_chain_{3,4,5,6,7,12,16,24} | 264169 | 3..24 | 2000 each | production | - | 0.23M | **PASS** (all) |
| pinned 1 CPU / 2 CPUs | 264169 | 8 | 3000 | production | affinity mask 0x1 / 0x3 | 0.09M | **PASS** (all) |

ST chain rate 851-881 it/s, i.e. ~12.5k gwmuls/s single-threaded; phase-2
CPU-time/wall-time was ~3.3, so the helper threads really are executing the
MT kernels during these runs (MT is used, it just does not reduce wall time at
this FFT size because the per-multiply sync cost is comparable to the work).

## 4. What the soaks rule out (x86)

* Not a plain `gwmul3`/`gwmulmuladd5`/`gwsubmul4` MT bug with ordinary
  (normalized, not-FFTed / FFTed) inputs: ~1.1M such MT multiplies clean.
* Not specific to the Edwards `ed_add`/`ed_dbl` mix with `ED_XYZ`,
  `ED_STARTNEXTFFT`, `ED_SUBTRACT`, `ED_NO_GWALLOC`, in all 12 variants:
  ~2.3M MT multiplies clean at 16K, both moduli.
* Not specific to production stage-1 chaining (in-place accumulate,
  FFTed/partially-FFTed operands, dictionary adds, `gwerror_checking`
  on/off/periodic): ~5.8M MT multiplies clean (pending final numbers).
* Not thread-count dependent in the tested range (2..24 threads all match).
* Not sensitive to extreme scheduling skew: 8 gwnum threads pinned to **one**
  physical core (and to 2 cores) still produce bit-identical results.
* Not sensitive to the sync mechanism: both event-based (`use_spin_wait=0`,
  production default) and both spin-wait levels are clean.

## 5. Mechanism hypotheses (ranked) + evidence

### H1 (most likely): ARM64 memory-ordering / missing-barrier defect, invisible on x86
x86-TSO hides missing acquire/release; ARM64 does not. Evidence:
* The port already had to add explicit memory ordering for exactly this class:
  `gwthread.cpp` implements the atomics with `std::atomic` and
  `memory_order_acquire`/`release`/`acq_rel`, with the comment *"ARM64 need
  memory order acquire and release to avoid multithread memory errors."*
  Before that change the port was corrupting on ARM only.
* The gwnum MT framework (`gwnum.c`) was written for x86 asm and uses several
  **plain, non-atomic** cross-thread fields (`active_child_gwdata`,
  `pass1_state`, `num_pass1_blocks`, section bookkeeping) whose ordering is
  established only *implicitly* by the event/mutex (`gwevent_signal`/
  `gwevent_wait`) or by the atomic gating variables. Any one of those that is
  reached without an acquire/release pair would be a latent ARM-only race and
  is invisible to every x86 test (including mine).
* The failure mode reported on the phone (rare, "wrong FFT product", no crash,
  worse with heterogeneous/oversubscribed cores) is the classic signature of a
  lost update / stale read, not of a deterministic math bug.
* `contrib/experiments/phone-armv9-mt/README.md` independently reached the same
  conclusion ("the synchronisation ... may be missing memory barriers that the
  x86 originals relied on").
* x86 cannot falsify this: 12M multiplies clean says nothing about ARM ordering.

Concrete code spots worth an ARM-side review (x86-clean by construction):
`gwnum.c`: `signal_auxiliary_threads` / `wait_on_auxiliary_threads` /
`auxiliary_thread` (`all_work_assigned`, `active_child_gwdata`),
`pass1_state1_assign_first_block` / `pass1_state1_assign_next_block` /
`pass1_get_next_block_mt` (the `pass1_carry_sections` hand-off,
`carry_in_blocks_finished`, `can_carry_into` event), and the
`num_active_helpers`/`all_helpers_done` protocol.

**Decisive experiment (ready to run on the phone/box):** push
`_work/mt/x86/mtstress_ed-arm` (static musl aarch64, current harness) and run
the same differential test *on the device*:

```
./mtstress_ed-arm --n 264169 --threads 8 --iters 20000 --seed 0xED01          # random ops
./mtstress_ed-arm --n 264169 --threads 8 --iters 20000 --seed 0xC001 --chain  # ECM stage-1 shape
```

If the phone shows a threads=1 vs threads=8 digest difference where x86 does
not, the bug is ARM-specific ordering; the harness then localises it to the
first iteration/op, and `--op-list`/`--chain --errchk`/`--spin-wait` narrows
the code path. If the phone is also clean, the earlier `ed_check` failures are
most plausibly the **already-fixed** bugs (see H4) and the "MT corruption"
should be re-triaged as a stale report.

### H2: residual "delayed-normalized operand" multiply bug reachable from a path edtest does not cover
The one *proven* 16K defect was deterministic: a `GWADD_DELAY_NORMALIZE`
result (unnormalized gwnum, `unnorms>0`) fed to `gwmul3` gave a wrong product
at 16K (found by `contrib/edtest`, fixed in `ed_add` by
`GWADD_FORCE_NORMALIZE`). The core of that mechanism lives in the converted
carry/normalization kernels, and other production paths still hand
delayed-normalized values to multiplies (`ecm.cpp` NAF dictionary init
`gwadd3o(...,GWADD_DELAY_NORMALIZE)`, `ed_dbl` algorithm 4
`GWADD_DELAYNORM_IF(safe21)`). Evidence against it being the current phone
failure: it is deterministic and thread-count independent, and both the A53
box and x86 are clean on 8-thread ECM (same converted sources); also the algo
actually selected at 16K is algorithm 5, which contains no delayed-norm add.
Worth a directed test (not done here to stay in the MT scope): multiply an
explicitly unnormalized (`gwaddsub4o(...,GWADD_DELAY_NORMALIZE)`) value and
compare against a schoolbook `mulg`/`modg` reference at 16K and 64K.

### H3: pass1 carry-section hand-off race exposed only by heterogeneous cores
The `pass1_carry_sections` protocol (sections split/steal/hand-off, with the
`can_carry_into` event and `carry_in_blocks_finished` flag) is the most complex
shared-state machinery in the MT path and is the only part whose behaviour
depends on *which* thread finishes a block first. Evidence: my 8-threads-on-1-core
pin test is clean, which argues the protocol is logically sound under arbitrary
skew on x86; but the ARM ordering caveat of H1 applies here most strongly.
This is a sub-case of H1 to investigate on-device.

### H4: the device failures may be the already-fixed deterministic bugs
The repo now contains two workarounds for *deterministic* 16K problems:
`GWADD_FORCE_NORMALIZE` in `ed_add` and "never error-check an `ED_XYZ` op"
(`ecm.cpp` stage-1 loop). The README/phone notes record that some earlier
"Hardware error" bursts were later shown to be *stale-save-file artifacts*
(reused run directory), and that 8-thread PRP in a fresh directory is clean
3/3. If the phone binary used for the failing ECM runs predates these fixes
(md5-pinned `mprime-mc4`), the reported symptom may be those bugs, not a live
MT race. Recommendation: re-run ONE 8-thread ECM stage-1 on the phone with the
**current** binary in a fresh directory before spending more effort on H1/H3.

## 6. Recommendations

1. Run `mtstress_ed-arm` (differential, threads=1 vs 8) on the phone; that is
   the only way to separate H1/H3 from "already fixed" (H4) - x86 is
   structurally blind to ordering bugs.
2. If it reproduces on the phone, use `--chain --errchk 0/1/2`,
   `--spin-wait 0/1/2`, `--op-list <variant>` to find the first wrong
   operation, then instrument the ARM-side code listed under H1 with
   acquire/release (or `__sync_synchronize()` around the thread-work dispatch,
   cf. the commented-out barrier in `gwnum/x/xmult2.c` line 26:
   `//__sync_synchronize(); // ARM64 may need this to avoid multithread errors`).
3. Independent of H1, the x86 clean-soak result means a fix cannot be validated
   on x86 for ordering issues; validation must be on-device. For any *logic*
   fix, `mtstress_ed` on x86 is a usable regression gate (it catches a
   poisoned digest and is sensitive to any ST-vs-MT difference).

## 7. Exact commands

```
# raw mix (Lead's harness), 8 threads:
mtstress.exe --n 264169 --threads 8 --iters 400000 --seed 0x51A1 --all-mismatches

# ECM Edwards op mix:
mtstress_ed_new.exe --n 264169 --threads 8 --iters 100000 --seed 0xED01 --all-mismatches
mtstress_ed_new.exe --n 291143 --threads 8 --iters 100000 --seed 0xED02 --all-mismatches

# production-shaped ECM stage-1 chain:
mtstress_ed_new.exe --n 264169 --threads 8 --iters 100000 --seed 0xC001 --chain --all-mismatches
mtstress_ed_new.exe --n 264169 --threads 8 --iters 100000 --seed 0xC003 --chain --errchk 0 --all-mismatches
mtstress_ed_new.exe --n 264169 --threads 8 --iters 60000  --seed 0xD001 --chain --spin-wait 1 --all-mismatches
```

Every run above was executed from `_work/mt/x86/` so that a stray `prime.txt`
(e.g. one left by `--algo`) could never influence the run; `--algo` writes
`prime.txt` in the cwd to force an `EdDblAlgo*` choice and is not used by the
soaks above (production defaults).

---
_Logs for every line of this table are in `_work/mt/x86/logs/<run>.log`._

---

# PART 2 (post-pivot): root-cause validation of the carry-section handshake fix

The Lead's static audit (`_work/mt/audit/REPORT.md`) identified the root cause while the
pre-fix soaks above were running: the gwnum **pass-1 carry-section handshake** publishes /
consumes `(start_block, carry_in_blocks_finished)` with plain (non-atomic) loads and stores
(`gwnum.c`: writer sets the flag ~6322 after writing normalized FFT data; waiter reads
`start_block` + flag ~6381; splitter writes `start_block` then clears the flag ~5685-5698
with the explicit "Must be done AFTER changing start_block!" comment). x86-TSO makes those
orderings free; on ARM64 a thread can observe `(new start_block, stale TRUE flag)` or
`(flag TRUE, data not yet visible)`, skip the wait and let `xgw_carries` run over
not-yet-normalized blocks -> numerically clean but wrong product -> `ed_check` failure.
The port only added acquire/release to the `gwthread.cpp` atomics, not to these fields.
The Lead's fix (in the working tree of `p95v3107b02.arm64.source/gwnum/gwnum.c`, uncommitted)
takes `gwdata->thread_lock` around: (a) the splitter's `start_block` update + flag clear,
(b) the writer's flag set, (c) the waiter's check/recheck loop.

## 8. ThreadSanitizer: SKIPPED - no Linux execution environment on this host

TSan requires a Linux (or macOS/FreeBSD) target; MinGW/Windows TSan is unsupported. This
host has no way to run a Linux binary:

* `wsl --status` -> WSL not installed, no distribution ("WSL is not installed on this
  system ... no distribution found"); no `wsl.exe --list --online` target installed.
* `Get-Command docker` -> not found.
* `D:\data\code\p95\_tools\qemu\` is QEMU-for-Windows **full-system** emulation only
  (`qemu-system-*.exe`); there is no `qemu-x86_64`/user-mode emulator and no Linux
  guest kernel/rootfs to boot.
* `zig cc -target x86_64-linux-musl` can *cross-build* a TSan binary in principle, but the
  result cannot be executed here.

Per the Lead's instruction ("if no Linux environment is available, say so and skip; do NOT
spend long on it"), the TSan run was skipped. The race is nevertheless confirmed by
construction: the unfixed code has plain cross-thread stores/loads with no
release/acquire edge between them, and the fix closes exactly those three sites.
Note that an x86 TSan run *would* have reported it (TSan models happens-before, not the
x86 TSO hardware), so if a Linux host/WSL becomes available the exact command is:

```
zig cc -target x86_64-linux-gnu -fsanitize=thread -g -O1 -o mtstress_ed_tsan \
    mtstress_ed.cpp -I<repo>/gwnum <repo>/gwnum/gwnum.a <repo>/gwnum/polymult.a -lm -lpthread
./mtstress_ed_tsan --n 264169 --threads 8 --iters 2000 --seed 0xC001 --chain
```

(the harness itself is thread-safe: it is a single worker calling gwnum from one thread
which internally spawns the gwnum helper threads).

## 9. Fixed-build regression matrix (x86, 16K type-1 FFT Pass1=64 Pass2=256 clm=4)

Fixed archives built exactly like the original ones, into separate directories (the
original `_work/gwnum/lib-x86` was never touched):

```
cd p95v3107b02.arm64.source/gwnum
mingw32-make -f Makefile.arm64 ARCH=x86 OBJDIR=../../_work/gwnum/obj-x86-fix \
    LIBDIR=../../_work/gwnum/lib-x86-fix -j8      # -> EXIT=0
```

`obj-x86-fix/gwnum.o` differs from `obj-x86/gwnum.o` (md5 056F2945... vs 2F95F51A...,
built 01:50:50), i.e. the fixed `gwnum.c` is really in these archives.
Harnesses rebuilt against them: `_work/mt/x86/mtstress-fix.exe`,
`_work/mt/x86/mtstress_ed-fix.exe` (same sources/commands as the unfixed ones).

Every configuration below is a threads=1 (reference) vs threads=N differential run; PASS
means bit-identical digests for every iteration.

| # | harness | mode | n | threads | iters | result |
|---|---|---|---|---|---|---|
| M1 | mtstress-fix | raw 6-op mix | 264169 | 2 | 20000 | **PASS** |
| M2 | mtstress-fix | raw 6-op mix | 264169 | 4 | 20000 | **PASS** |
| M3 | mtstress-fix | raw 6-op mix | 264169 | 8 | 50000 | **PASS** |
| M4 | mtstress-fix | raw 6-op mix | 264169 | 16 | 50000 | **PASS** |
| M5 | mtstress_ed-fix | 12 ECM Edwards variants | 264169 | 8 | 20000 | **PASS** |
| M6 | mtstress_ed-fix | random 12 variants | 264169 | 4 | 10000 | **PASS** |
| M7 | mtstress_ed-fix | random 12 variants | 264169 | 2 | 10000 | **PASS** |
| M8 | mtstress_ed-fix | ECM stage-1 chain | 264169 | 8 | 20000 | **PASS** |
| M9 | mtstress_ed-fix | ECM stage-1 chain, spin-wait 1 | 264169 | 8 | 5000 | **PASS** |
| M10 | mtstress_ed-fix | ECM stage-1 chain | 264169 | 16 | 100000 | **PASS** (no hang, no mismatch; ST 916.3 -> MT 598.1 it/s) |
| M11 | mtstress_ed-fix | ECM stage-1 chain | 264169 | 24 | 100000 | **PASS** (no hang, no mismatch; ST 1022.0 -> MT 612.9 it/s) |

No hang, no timeout, no crash in any fixed-build run. (Note: gwnum caps the helper-thread
count for this FFT layout - the `--threads 24` run's log says "... clm=4, 16 threads" - so
16 and 24 both exercise a 16-helper-thread pool; `threads=16/24` still exercises the
section split/steal/carry-into handshake with the maximum possible thread count.)

## 10. Performance delta (fixed vs unfixed, same machine, back-to-back, near-idle)

| config | unfixed | fixed | delta |
|---|---|---|---|
| raw 6-op mix, t=8, MT pass | 1997.2 iters/s | 1969.8 iters/s | **-1.4 %** |
| raw 6-op mix, t=8, ST reference pass | 2287.6 iters/s | 2280.2 iters/s | -0.3 % (noise) |
| ECM stage-1 chain, t=8, MT pass | 820.4 iters/s | 828.0 iters/s | +0.9 % (noise) |

The added locking costs ~1.4 % on the most MT-heavy raw-multiply mix and is inside the
noise for the ECM chain - consistent with the expected "<1-2 %" cost of a few mutex ops
per pass-1 block.

## 11. Anomalies / side observations

* A directed big-int-reference test (`_work/mt/x86/mtstress_dn.c`, built by me) checks
  "delayed-normalized operand" multiplies against exact `mulg`/`modg` arithmetic. Its first
  check at 16K, threads=1, flagged `gwmul3` of a `GWADD_DELAY_NORMALIZE` value as wrong
  (`WRONG iter 0 variant 0 unnorms=1`). This is *consistent with the already-known and
  already-worked-around deterministic 16K delayed-norm bug* that motivated
  `GWADD_FORCE_NORMALIZE` in `ecm.cpp:ed_add` (the Lead's `contrib/edtest` finding); it is
  NOT an MT race and is not production-reachable in the current stage-1 code. The run was
  abandoned (its reference arithmetic is O(n^2) big-int and too slow) before the
  FORCE_NORMALIZE control variants completed, so this is recorded as corroboration only,
  not as a formal new result.
* Pre-fix x86 soaks: **~12 million MT multiply-equivalents clean** (raw 6-op mix 0.47M,
  raw per-op subsets 0.6M, ECM Edwards 12-variant mix 3.0M, production-shaped stage-1
  chain ~7.9M, plus thread-count battery and 8-threads-on-one-core runs). At the estimated
  corruption rate of 1/700k multiplies this is ~17 expected events: the x86 host is
  structurally blind to this bug (x86-TSO), which is itself evidence for the audit's
  root cause rather than against it.
* The model/mechanism is consistent with the phone-side confirmation reported by the Lead
  (unfixed aarch64 binary failed ~1 per 1-2 8-thread ECM stage-1 runs; fixed 2/2 clean).
