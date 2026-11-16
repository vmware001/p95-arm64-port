# contrib/experiments/mt-fix — the ARM64 `CoresPerWorker > 1` corruption bug (fixed)

**Status: FIXED.** Per-worker multithreading on ARM64 intermittently returned wrong
FFT products. Reproducer: one worker, 8 compute threads (`NumCPUs=8`, `NumWorkers=1`,
`CoresPerWorker=8`) running ECM stage 1 on `M264169`, `B1=50000` (16K FFT,
`Pass1=64, Pass2=256, clm=4`); stage 1 reaches 100.00 %, then `ed_check` fails:

```
[Worker] Hardware error?  Point no longer on the Edwards curve.  Restarting from last save file.
```

Pre-fix failure rate on a Dimensity 9400 phone: roughly **one error per 1–2 full
stage-1 runs** (72,121 exponent bits per run, ~145 s per run at 8 threads). The same
signature was reported from Termux on the same phone (45 %/80 % of stage 1 depending
on B1) and is the reason `CoresPerWorker=1` was previously mandatory.

* Single-threaded runs were never affected (the farm has run 4×1 threads for days,
  and 12 h+ per ECM curve on the RK3528 box).
* 8-thread **PRP** (squares only) looked clean; ECM does ≈10 multiply-class FFT
  passes per ladder bit vs 1 square per PRP bit, so it hits the bug ~10× more often.
* 4-thread ECM looked clean in ad-hoc tests; with 8 threads gwnum has at most ~5
  carry sections (`num_pass1_blocks = 128`, `min_split_size = 24`), so helpers run
  out of work and **steal/split other threads' live sections**, which is the code
  path that races.

## Root cause (ARM64 memory model, not the converted kernels)

`gwnum/gwnum.c`'s pass-1 carry-section state machine coordinates the compute
threads with **plain loads/stores** on `struct pass1_carry_sections`
(`gwnum.h:923-931`): `carry_in_blocks_finished`, `start_block`, `last_block`,
`next_block`, `carry_out_section`, `section_state`. Upstream gwnum is correct on
x86-64 because x86-TSO makes stores from one core visible to other cores in program
order. The source even states the requirement:

* `gwnum.c:5689` — `carry_in_blocks_finished = FALSE; // Must be done AFTER changing start_block!`
* `gwnum.c:6344-6345` — "SECTION SPLITTING CODE MUST CHANGE THESE VALUES IN A
  PARTICULAR ORDER SUCH THAT THIS CHECK ALWAYS WORKS."

ARM64 has a weak memory model, so two things can happen on a helper thread:

1. the waiter can observe a section's **new**
   `start_block` together with a **stale** `carry_in_blocks_finished == TRUE`
   left over from that section slot's previous use (splitting re-uses slots), skip
   the `can_carry_into` wait, and run `xgw_carries` over the section's first blocks
   *before* the owning thread has finished normalising them; or
2. see `carry_in_blocks_finished == TRUE` before the preceding normalisation data
   writes are visible (the flag store has no release semantics; the event signal is
   only reached on the slow path).

Either way the carry is added to pre-normalisation inverse-FFT values. The
following forward FFT then spreads a sub-ULP difference over the spectrum: the
product is **wrong but numerically clean** (roundoff checking never fires), so the
error only surfaces much later at the next `ed_check` — exactly the observed
"stage 1 finished, then restart from the last save file".

The converted FFT/multiply/normalisation kernels themselves were audited and are
clean (see `REPORT-static-audit.md`: no static/shared state, no `thread_num`
dependence, faithful asm translation of the carry/normalisation macros, MT section
indexing confined to `gwnum.c`).

## The fix (`gwnum/gwnum.c`, commit "gwnum: fix ARM64 MT carry-section race")

Serialise the `(start_block, carry_in_blocks_finished)` handshake on gwnum's
existing `thread_lock` mutex, so the pair is coherent and the normalised data is
published before the flag:

1. **section splitter** (`pass1_state1_assign_first_block`): the new section's
   `start_block/next_block/last_block/carry_out_section` writes and the
   `carry_in_blocks_finished = FALSE` clear now happen under `thread_lock`;
2. **flag setter** (`pass1_get_next_block_mt`, end of the first-blocks
   normalisation): sets `carry_in_blocks_finished = TRUE` under `thread_lock`
   (mutex release publishes the FFT data written just before it);
3. **carry-into waiter** (`pass1_get_next_block_mt`): reads the pair under
   `thread_lock`, keeping upstream's counter/event lost-wakeup dance intact.

These paths only execute when `num_threads > 1` (the single-thread path uses
`pass1_get_next_block`, not the `_mt` variant), so single-threaded farms are
bit-for-bit unaffected. Cost: a few uncontended mutex operations per pass-1 block
per thread — below measurement noise.

## Verification

| test | result |
|---|---|
| phone, **unfixed** mprime-mc4 (`f630c8ac…`), 8-thread ECM M264169 B1=50000, full stage 1 | **FAIL** — `Point no longer on the Edwards curve` at 100.00 %, restart from save (`phone-repro-unfixed-excerpt.txt`) |
| phone, **unfixed** mprime-mc4, 8-thread ECM M291143 B1=50000 (16K FFT, `3.019 bits-per-word below FFT limit` = the owner's original Termux case) | **FAIL** — 1 hardware error in 1 run (`summary-unfixed-291143.txt`); this is the exact signature from the screenshot |
| phone, **fixed** binary, 8-thread ECM M264169 B1=50000, full stage 1, **5 runs** | **5/5 PASS** — reached stage 2, 0 roundoff/hardware errors (145/145/140/145/145 s per run) |
| phone, **fixed** binary, 8-thread ECM M291143 B1=50000, **2 runs** | **2/2 PASS** — reached stage 2, 0 errors (140/145 s) |
| x86 differential MT-vs-ST stress harness (`mtstress.c`, 1e6 iters, 8 threads) | PASS before and after (x86-TSO hides the race; useful as a regression gate) |
| x86 fixed-build regression suite (11 configs: raw 6-op mix, ECM-variant `ed_add/ed_dbl` mix, production-shaped stage-1 chain; threads 1/2/4/8/16/24; spin-wait on/off) | **all PASS, no hangs**; 100k-iteration 16/24-thread chain soaks clean; perf −1.4 %…+0.9 % (noise, <2 % expected mutex cost) |
| x86 **unfixed** pre-fix soak (~12 M MT multiply-equivalents: raw/op-subset/ECM-mix/chain, threads 2–24, 8 threads pinned to one physical core) | **clean** — ~17 corruption events were statistically expected at 1/700k if x86 could see the race; this is itself evidence that the bug is the ARM64 memory-ordering hole (x86-TSO hides it). ThreadSanitizer was not runnable on the Windows build host (no WSL/Docker/qemu-user); the exact recipe is in `_work/mt/x86/REPORT.md` §8. |
| RK3528 box (4× A53), **unfixed** binary: 8-thread ECM stage 1 and the MT stress harness at 8/16 threads | did not reproduce in the available time (the box was never previously run to a full 8-thread stage 1; its earlier "clean" claim was only 8 % of a stage — see note below). The box is homogeneous A53, i.e. far less thread-speed skew than the phone, so the race window is much smaller there. |
| RK3528 box, **fixed** binary, 8-thread ECM M264169 B1=50000 | **PASS** — full curve completed end-to-end: *Stage 1 complete. 1,221,525 transforms … 560.9 s*, *Stage 2 complete. 31,910 transforms … 17.1 s*, `M264169 completed 1 ECM Edwards curve`, 0 Edwards/hardware errors (`box-fixed-8thread-curve-excerpt.txt`); this also retires the old "8-thread ECM hangs at polymult helper threads" observation |
| RK3528 box, **fixed** binary, LL known-answer harness `gwlltest --threads 8 --extended` | **ALL PASS** — 20 Lucas-Lehmer tests: every known Mersenne prime (127…216091, including 8K/12K FFTs at 8 threads) gives residue 0 and 2047 gives nonzero (`box-gwlltest-8thread-extended.log`) |
| RK3528 box, fixed binary, MT differential stress harness at 8/16 threads | PASS (also PASSed with the unfixed binary — like on x86 the MT harness is not sensitive enough to reproduce this <2e-6/transform race; it is kept as a regression gate) |
| static audit of converted kernels (`REPORT-static-audit.md`) | no kernel-side MT defect found; root cause is the `gwnum.c` handshake |

Notes:

* A full 8-thread stage-1 run is 72,121 exponent bits ≈ 145 s on the Dimensity
  9400 and ≈ 520 s on the RK3528 box (the box's earlier "8-thread ECM clean"
  result in `contrib/experiments/phone-armv9-mt/README.md` only reached 8.07 % of
  stage 1 because the farm was still running — it was not evidence either way).
* `mtstress.c` did not reproduce the corruption on the phone either (20 k and
  30 k iterations = 30–50 k multiplies, ~20× less exposure than one ECM
  stage 1); the corruption needs the ECM-style multiply mix *and* the section
  stealing that happens with 8 threads. The mprime ECM reproducer above is the
  reliable one.


## Residual risk / follow-up

The fix covers the `(start_block, carry_in_blocks_finished)` publication that the
race analysis identified. Other fields of the same state machine
(`section_state`, `next_block`, `last_block`, `carry_out_section`) are still
plain `int`s and are read by other threads without the change-in-progress lock
(see S2 in `REPORT-static-audit.md`); upstream explicitly tolerates that layout
on x86, and 7/7 verified phone runs plus the x86 regression suite show no
remaining corruption. If any MT anomaly ever reappears, apply the same
discipline (writes under `thread_lock`/release, reads after acquire) to those
fields, starting with the splitter's unlocked scan in
`pass1_state1_assign_first_block` and the `section_state` transitions in
`pass1_get_next_block_mt`.

Separately, `gw_get_maxerr()` is still neutralised on ARM64 (see `f195b22`) -
that is a *reported roundoff value* workaround, not an arithmetic change, and is
unrelated to this bug.

## Files

* `mtstress.c` — differential multithreaded stress harness (public gwnum API): runs
  a deterministic sequence of `gwmul3`/`gwmulmuladd5`/`gwsubmul4`/chained
  `STARTNEXTFFT` operations twice (threads=1 reference vs threads=N) and compares a
  digest of every result; prints the first mismatching iteration. It forces the x86
  build to use the same "type-1 FFT 16K, Pass1=64, Pass2=256, clm=4" layout and CPU
  flags as the ARM64 build. Build lines are in the header comment.
* `REPORT-static-audit.md` — the static audit that identified the `gwnum.c`
  carry-section handshake as the root cause (and ruled out the converted kernels).
* `phone-repro-unfixed-excerpt.txt` — excerpt of the failing phone run with the old
  binary (worker startup, stage-1 100 % + ed_check failure + restart).
* `phone-runfix.sh`, `phone-runfix2.sh` — the phone test scripts used to reproduce
  and to validate the fix (fresh directory per run, detached, cleaned up afterwards).
