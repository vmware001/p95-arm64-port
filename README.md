# p95-arm64-port

Prime95 / mprime **v31.7 build 2** ported to Linux **ARM64 (AArch64)** — validated on a
Rockchip RK3528 set-top box (4× Cortex-A53 @2.0 GHz, 2 GB RAM, Android 9 + chroot-ish
userland). See [README-arm64.md](README-arm64.md) for the full port notes.

> Upstream source is © Mersenne Research, Inc. / George Woltman and is publicly
> distributed at https://www.mersenne.org/download/ (`p95v3107b02.source.zip`).
> This repo keeps that code intact and adds the ARM64 port layer on top — same license
> terms as upstream apply to the upstream files.
> The ARM64 gwnum base is derived from [sam8641/gwnum_c1](https://github.com/sam8641/gwnum_c1)
> (SSE2→NEON vector-C conversion).

## What this repo adds

* `gwnum/Makefile.arm64`, `gwnum/arm64_sse2.h` — ARM64 build of gwnum (zig cc,
  `aarch64-linux-musl`, static)
* `linux-arm64/` — console mprime for ARM64 Linux
* **Bug fixes found on real hardware** (all with root-cause analysis in commit messages):
  * `72bb02a` ECM `ed_add` forced normalization of E/H — the NEON-converted gwnum returns
    a wrong product for the delayed-norm multiply at ≥16K FFT (proved by the independent
    checker below; this was the "30-minute restart wall")
  * `f195b22` neutralized the spurious MAXERR=0.5 readout on ARM64
  * `6fd0572`/`72fe6b2`-line fixes: P-1 stage-2 polymult layout, ECM stage-2 SIGSEGV,
    PrimeNet `ss`/`sh` parameters
  * **per-worker multithreading** (`CoresPerWorker` > 1) — gwnum's pass-1 carry-section
    handshake relied on x86-TSO store ordering; on ARM64 a thread could propagate
    carries into not-yet-normalized FFT blocks and produce a wrong (but numerically
    clean) product, which showed up as `Point no longer on the Edwards curve` after a
    full ECM stage 1.  Fixed by serialising the `(start_block, carry flag)` handshake
    under gwnum's `thread_lock` — see `contrib/experiments/mt-fix/` for the analysis,
    reproducer and verification
* `contrib/primenet/` — **manual PrimeNet workflow** scripts (source-built clients are
  excluded from the automatic v5 protocol: `pnErrorResult=9`; these drive the web pages
  instead — login / assign / submit + an auto-farm loop for the RK3528 box)
* `contrib/edtest/` — independent Edwards-arithmetic checker (links only gwnum, verified
  against Python big-int math) that pinned the 16K-FFT `ed_add` bug
* `contrib/experiments/mt-fix/` — the **`CoresPerWorker>1` corruption bug**: static
  audit, phone reproducer, differential MT stress harness (`mtstress.c`) and the
  before/after validation table

## Build

Windows host, zig 0.15.1 cross-compile:

```
# gwnum library
cd gwnum && mingw32-make -f Makefile.arm64 ZIG=/path/to/zig

# mprime binary (see _work build logs in the repo history)
cd linux-arm64 && zig c++ -target aarch64-linux-musl -static -s -o mprime \
    prime.o menu.o cJSON.o ecm.o exponentiate.o pair.o pm1prob.o factor64_c.o \
    ../gwnum/gwnum.a ../gwnum/polymult.a \
    libhwloc.a libcurl.a libgmp.a -lm -lpthread
```

## Hardware status

| Work type | Status |
|---|---|
| LL / PRP / TF | validated correct (14/14 LL residues, Gerbicz OK) |
| ECM | **fixed** — was broken at ≥16K FFT by the delayed-norm bug; now curve-complete |
| P-1 | stage 1+2 OK (tiny exponents have a known save-path quirk) |

PrimeNet: source builds are *untrusted* (server policy) — use `contrib/primenet/` for
assignment/result automation, or run offline.

## Prebuilt binary

Grab the static AArch64 build from
[Releases](https://github.com/vmware001/p95-arm64-port/releases):
**`mprime-mc4mt`** (10,126,600 bytes, md5 `c82f4fbcf27c2ac0cbb71ba7f8337088`) —
includes every fix below, in particular the `CoresPerWorker>1` multithreading fix.
It is musl-static — runs on any arm64 Android (Termux, no root needed) or arm64
Linux with zero dependencies. (`mprime-mc4` / `mprime-mc3` are the previous
builds without the MT fix.)

### Termux quick start (no root)

```sh
pkg update && pkg install wget python
mkdir -p ~/mprime && cd ~/mprime
wget https://github.com/vmware001/p95-arm64-port/releases/download/v3.10.7b02-arm64/mprime-mc4mt -O mprime
chmod +x mprime

# minimal config BEFORE first run (in-client PrimeNet is unusable on untrusted builds:
# the v5 API answers pnErrorResult=9, so keep UsePrimenet=0 and drive the web flow manually)
cat > prime.txt <<'EOF'
NumCPUs=8
NumWorkers=8
Memory=2048
ProofPower=6
UsePrimenet=0
EOF

./mprime            # starts computing worktodo.txt; Ctrl+C exits, progress is saved
termux-wake-lock    # keep CPU awake
```

Work management (assign / submit) with the bundled manual client
(`contrib/primenet/primenet_manual.py`):

```sh
wget https://raw.githubusercontent.com/vmware001/p95-arm64-port/master/contrib/primenet/primenet_manual.py
export PN_PASS='your-primenet-password'   # never hardcode; script never writes it to disk
python primenet_manual.py --user YOURUSER login
python primenet_manual.py --user YOURUSER assign --pref 5 --exp-lo 100000 --exp-hi 400000 >> worktodo.txt
# ... after curves finish (results land in results.json.txt):
python primenet_manual.py --user YOURUSER submit --file results.json.txt
```

Termux survival tips: run inside `~/` (not `/sdcard` — noexec + slow), exempt
Termux from battery optimization in Android settings, and expect phones to
thermally throttle (a TV box/stick with a heatsink sustains full clocks far
better). ECM B1 must match the assigned bounds exactly or the server silently
drops the result.

### Fixed: per-worker multithreading (`CoresPerWorker` > 1)

Per-worker multithreading used to intermittently return wrong FFT products on
ARM64. A single ECM worker using 8 threads on a 16K FFT (e.g. `M264169`,
`B1=50000`) would finish stage 1 and then trip
`Hardware error? Point no longer on the Edwards curve. Restarting from last save file.`
about once every 1–2 stage-1 runs. The single-threaded path was never affected.

Root cause: `gwnum/gwnum.c`'s pass-1 carry-section state machine coordinates the
helper threads with plain cross-thread fields and relies on x86-TSO store
ordering (the source itself says *"SECTION SPLITTING CODE MUST CHANGE THESE
VALUES IN A PARTICULAR ORDER SUCH THAT THIS CHECK ALWAYS WORKS"*). On ARM64 a
thread can observe a section's new `start_block` together with a stale
`carry_in_blocks_finished = TRUE` from that slot's previous use, skip the wait,
and add carries into FFT blocks that have not been normalized yet. The result is
a wrong but numerically clean product — round-off checking never fires — so the
corruption only surfaces at the next `ed_check`.

Fix: serialise the `(start_block, carry_in_blocks_finished)` handshake on
gwnum's existing `thread_lock` mutex and publish the normalized FFT data before
setting the flag. Single-threaded mode is untouched (these paths only run when
`num_threads > 1`). Validation: 5/5 clean full 8-thread ECM stage-1 runs on the
fixed binary, the unfixed binary fails the owner's original `M291143` case
(`3.019 bits-per-word below FFT limit`, 8 threads) in 1 run, plus the box and
LL-harness regression tests — full analysis, reproducer and numbers in
[`contrib/experiments/mt-fix/`](contrib/experiments/mt-fix/README.md).
