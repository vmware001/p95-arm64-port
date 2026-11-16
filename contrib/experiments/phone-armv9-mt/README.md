# contrib/experiments/phone-armv9-mt — multithreaded FFT on a heterogeneous ARM phone

**Status: single-thread is proven exact and fast; multithreaded PRP completes
cleanly in a fresh directory (3/3 runs, 8 threads, full tests). An earlier
264,177-line "Hardware error" burst turned out to be a *test-harness artifact*
(reusing a work directory that still contained an ECM save file), not a real
MT defect. 8-thread ECM still hangs at polymult helper-thread setup - that is
the remaining open question. Until it is settled, production phones should run
`CoresPerWorker=1`.**

## Device under test

| | |
|---|---|
| model | OPD2506 (OnePlus/OPPO), SoC `mt6991` = MediaTek Dimensity 9400 |
| CPU | 8 cores, heterogeneous: 1x Cortex-X925 + 3x Cortex-X4 + 4x Cortex-A725 (ARMv9) |
| OS | Android 17 (SDK 37), kernel 6.6.118-android15 |
| access | wireless adb (`adb pair` + `adb connect`), `u:r:shell:s0`, **no root** |
| binary | the same validated `mprime-mc4` used on the box, md5 `f630c8ac00088d36aa8e54e7e4811500` |

The box for comparison: RK3528, 4x Cortex-A53 @2016 MHz, homogeneous.

## Measured results

| test | result |
|---|---|
| PRP M11213 + M1277, 1 thread pinned to CPU7 | `probable prime` / `not prime` correct, Gerbicz 4/4 pass, 0 errors |
| **ECM M264169 B1=50000, 1 thread pinned to CPU7** | stage 1 **complete**: 72,126 exponent bits in **~108 s**, 0 errors, then stage 2 (`125 MB, D=1050`) |
| **PRP M264169, 8 threads (all cores), fresh dir, run to completion** | **clean, 3/3 runs**: completes in ~60 s, `Hardware error` x0, `Gerbicz error check passed` x2 |
| same, restricted to the A725 cluster (mask `0f`) | clean, completes in ~60 s |
| PRP 8 threads in a directory that still held an ECM save file | **`Hardware error` x264,177** and `ERROR: Comparing Gerbicz checksum values failed. Rolling back to iteration 0.` — **artifact of the reused directory** (see "Corrections" below) |
| **ECM M264169, 8 threads (all cores)** | **hangs** at `Setting affinity to run polymult helper thread on CPU core #8` (no progress for 120 s, killed) - still open |
| 2 workers x 1 thread, ECM M264169 + M279863 | stage 1 clean and fast (M264169: 72,126 bits in ~130 s with two workers running); stage 2 allocates 498 MB per worker and the phone starts swapping (3.9 GB swap used), so stage 2 progress is silent and slow |
| box, 8-thread PRP M264169, 16K FFT | clean (`Gerbicz error check passed`, no hardware errors) |
| box, 8-thread ECM M264169 B1=50000 | clean (4 min, 5817 stage-1 iterations, no Edwards/hardware errors) |

The Termux failure reported earlier (`Hardware error? Point no longer on
Edwards curve`, ~80% of stage 1 with an 8-thread worker) is the same defect:
the corrupted FFT makes the elliptic-curve point coordinates invalid.

## Stage 2 is silent, not hung (important)

ECM stage 2 logs almost nothing between "Using XXX MB of memory ... Ftree polys
in memory" and the curve completion - a run can look frozen for many minutes.
Sampling `/proc/<pid>/stat` (utime+stime) while a 4-thread ECM test sat on that
line showed the process CPU time **rising steadily for 6+ minutes** (no stall),
so it is computing, not deadlocked. A watchdog that only looks at log growth
will wrongly report a hang. Use CPU-time sampling and wait for
`completed 1 ECM curve` / `has a factor`.

Status of multithreaded ECM on the phone (fresh directory, 16K FFT):

| threads | behaviour |
|---|---|
| 4 | stage 1 completes (~100-130 s), stage 2 starts (polymult helper threads on cores 2-4), long silent CPU-bound stage 2 - no errors, no deadlock observed |
| 8 | stage 1 completes, then the curve **restarts from 0 %** and logs 2 `Hardware error` lines - suspicious, needs a longer fresh retest before use |

Practical guidance: for farm throughput, several single-threaded workers
(`CoresPerWorker=1`) are both safer and faster than one multithreaded worker
(MT gives ~2-3x on 4 threads, while 4 single-thread workers give ~4x);
multithreading is only attractive for using cores beyond what single-thread
workers already saturate.

## Corrections (important - two false conclusions were made and fixed)

1. **The 264,177 "Hardware error" PRP run was an artifact.** The four-phase
   test driver ran all phases in one directory and only deleted `*.bu` backups,
   but ECM save files are named `e<exp-hex>` *without* an extension. The PRP
   phase for the same exponent then found that stale file and produced garbage,
   which mprime reported as "Hardware errors!" plus a Gerbicz rollback loop.
   Re-running the same 8-thread PRP in a **fresh** directory: clean 3/3.
   Lesson: one fresh directory per test/phases, and match save files by the
   `e*` prefix, not by extension.
2. **The `Edwards curve` failures seen earlier in Termux were not analysed with
   a clean directory / current binary either** - they need a fresh reproduction
   before being blamed on multithreading.

## What this points at

1. **Heterogeneous cores / scheduling.** gwnum's multithreaded pass-2 relies on
   helper threads staying in lockstep (spin-wait barriers). With a 2.3 GHz
   X925 next to 2.0 GHz A725s the helper threads drift; the synchronisation in
   the ARM64 "asm converted to C" port (gwnum_c1) may be missing memory
   barriers that the x86 originals relied on (strong x86 memory model) -
   homogeneous A53 timing apparently hides it.
2. **mprime core numbering.** On this device mprime tries to pin a helper
   thread to "CPU core #8" although the cores are numbered 0-7, so thread
   placement can be wrong or fail (visible in the ECM hang).

Homogeneous-cluster affinity experiments (restricting the whole process to the
A725 cluster, the X4 cluster, etc.) are in `mt.sh` and are the next step: if
multithreading becomes clean inside one cluster, the trigger is heterogeneity
and the practical workaround is cluster pinning; a real fix would be adding the
missing barriers in the MT synchronisation code.

## Practical configuration for a heterogeneous phone today

```
NumCPUs=1-4          (one core per worker)
NumWorkers=1-4
CoresPerWorker=1     <- mandatory; any value > 1 hits the broken MT path
Memory=384-512       (shared by the workers; ECM stage 2 on a phone is
                      memory-hungry and starts swapping)
UsePrimenet=0        (manual submission only, see contrib/primenet)
```

Single-thread performance is excellent: an ECM B1=50000 stage 1 that takes the
A53 box ~50 minutes (under 4-worker load) completes in ~108 s on one X925 core.

## Scripts in this directory

- `mt.sh` — multithread stress test: `sh mt.sh <threads> <cpu-mask> <secs>`,
  where the mask is a taskset bitmask (Android/toybox `taskset` has no `-c`;
  use e.g. `ff` all cores, `0f` cpu0-3, `70` cpu4-6, `f0` cpu4-7).
- `run.sh` — the single-thread/8-thread four-phase validation used on the phone
  (PRP known answers, ECM single-thread, PRP 8 threads, ECM 8 threads).
- `ec2.sh`, `ec2b.sh` — two single-threaded ECM workers (the farm's real
  shape), with the "fresh directory" fix for stale ECM save files.
- `box_ecm8.sh` — the same 8-thread ECM test run on the RK3528 box.
- `phone-val` — one-line Termux installer/runner (downloads the binary, embeds
  all configs, runs the four phases, writes `~/phone-val/report.txt`).

Gotchas learned the hard way:
- ECM save files are named `e<exp-hex>` with **no** extension, so
  `rm -f ./*.bu` does not delete them (they are only periodic backups); use a
  fresh directory per run.
- Two identical worktodo lines are de-duplicated, so a second worker idles.
- `taskset -c N` does not exist in Android toybox; it is a bitmask.
- `timeout` sends SIGTERM and mprime's graceful shutdown can hang on Android
  17; use sleep + `kill -9` in test drivers.
