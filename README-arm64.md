# Prime95 / mprime for Linux ARM64 (aarch64)

This directory is the **ARM64 port** of Prime95 (mprime), built from the official
**Prime95 v31.07b02** source tree (gwnum **31.7**).

The original x86-only tree (`../p95v3019b20.source`, Prime95 v30.19b20) cannot be built for
ARM64: the gwnum FFT library is hand-written x86-64 assembly and no ARM64 version of gwnum 30.19
exists.  For gwnum 31.7 there is a community port of that assembly to portable C/vector intrinsics
([`sam8641/gwnum_c1`](https://github.com/sam8641/gwnum_c1), gwnum 31.7), which is used here.  That is
why this tree is based on the version-matched upstream **v31.07b02** sources instead of v30.19b20.

Everything below has been built and executed on a real (emulated) aarch64 Linux: see
"Validation" – the binary is a **static musl aarch64 executable** that runs LL, PRP, P-1 and
trial-factoring work.

---

## Status

| Component | State |
|---|---|
| mprime (CLI, aarch64) | **works** – static ELF64 AArch64, no interpreter, runs LL / PRP / P-1 / TF / benchmark |
| gwnum 31.7 FFT library | **works** – converted asm→C, 128-bit SIMD (NEON) kernels only |
| LL tests | **validated** – 12/12 known Mersenne primes give residue 0, composites correctly reject |
| PRP tests (Gerbicz) | **validated** – M2203 probable prime, M67 not prime; Gerbicz checks pass |
| Trial factoring (TF) | **validated** – C port of `factor64.asm`; factor found / no-factor cases match reference |
| P-1 factoring | **runs** (stage 1+2); stage 2 uses scalar C polymult – correctness not independently validated |
| ECM | links and runs, but **untested** (uses the same scalar polymult stage 2) |
| PrimeNet | links against HTTP-only libcurl (PrimeNet uses `http://`, matching upstream `--without-ssl`); not exercised from this machine |

### Known limitations

* **Only 128-bit SIMD kernels exist** on ARM64 (gwnum's SSE2/SSE4.1-class paths mapped to NEON).
  There is no AVX/AVX2/FMA3/AVX-512 equivalent, so large-FFT throughput per clock is below a modern
  x86-64 with the widest vector units.  No performance numbers from real ARM64 hardware are claimed –
  the only ARM64 machine used here was QEMU (TCG), which is *not* representative of speed.
* **P-1 / ECM stage 2 numerics are unvalidated**: only the scalar ("dbl", upstream 32-bit) polymult
  implementation exists.  T1 made it link- and layout-consistent (scalar dispatch, masked private CPU
  flags), but it has not been checked against the 128-bit gwnum data layout with a reference test.
  Use P-1/ECM results with that in mind; LL/PRP/TF are unaffected.
* **Trial factoring is ~10–100× slower** than the x86-64 assembly version (correctness over speed;
  GIMPS TF is normally run on GPUs anyway).
* libcurl is **HTTP-only** (no TLS backend).  This matches upstream's own `--without-ssl` builds and
  PrimeNet v5 URLs are plain `http://`; an `https://` URL would fail.
* hwloc is a minimal build (topology/caches/NUMA/bitmaps/binding from sysfs; no io/pci/GPU/libnuma).
* GMP is built with generic C `mpn` (no ARM64 assembly) – fine, GMP is not on the hot path.
* The security/CRC fields in results lines (`Wh8: 00000000,...`) can be all-zero on ARM64 because the
  CPUID/CPU-signature-derived values are synthetic; PrimeNet submission from ARM64 has not been tested.
* Windows ARM64 / macOS ARM64 (Apple Silicon) builds are **not** in scope.

---

## Building

### 1. Native build on an aarch64 Linux machine (recommended for users)

Requirements (Debian/Ubuntu names): `gcc g++ make libgmp-dev libhwloc-dev libcurl4-openssl-dev libboost-dev`
(RHEL/Fedora: `gcc gcc-c++ make gmp-devel hwloc-devel libcurl-devel boost-devel`).
Boost headers are required by `pair.cpp` (`<boost/circular_buffer.hpp>`), header-only.

```sh
cd gwnum
make -f Makefile.arm64 -j$(nproc)          # builds gwnum.a + polymult.a for aarch64
cd ../linux-arm64
make -f Makefile.arm64 -j$(nproc)          # builds ./mprime
./mprime -v
```

### 2. Static cross build with zig (used to produce the shipped self-contained binary)

This is what was used on a Windows x86-64 host to generate the tested artifact; it also works from
Linux/macOS hosts that have `zig` and `make`.

```sh
# 1. cross-built dependency archives (gmp, hwloc, curl) for aarch64-linux-musl
#    expected in ../../_deps/arm64/{lib,include}  (see ../PORTING-PLAN.md / _deps/README-deps.txt)
# 2. build
cd gwnum        && make -f Makefile.arm64 ZIG=/path/to/zig
cd ../linux-arm64 && make -f Makefile.arm64 CROSS=1 STATIC=1 ZIG=/path/to/zig \
     DEPS=../../_deps/arm64 EXE=../../_work/mprime/mprime-aarch64
```

Convenience wrappers: `linux-arm64/build-arm64-static.sh` (POSIX shell) and
`linux-arm64/build-arm64-static.ps1` (Windows PowerShell).  Both are relative-path based and
self-verify that the result is an ELF64 AArch64 executable.

The static binary is musl-linked and fully self-contained (no dynamic loader, no runtime deps).

---

## Running

```sh
./mprime -v                # version
./mprime -m                # interactive menu / first-time setup
./mprime -d                # daemon (no menu)
./mprime -t                # torture/stress test
./mprime -b                # benchmark
```

Typical first-run configuration (also used for the automated tests) – `prime.txt`:

```
UsePrimenet=0          # offline; omit to join PrimeNet
StressTester=0         # suppress the first-run wizard when running unattended
SkipTrialFactoring=1   # avoids the pre-LL trial-factoring pass (needed for small exponents)
NumCPUs=1              # optional: worker/thread override; default comes from hwloc
```

`worktodo.txt` examples (all validated on aarch64):

```
Test=127,0,1              # Lucas-Lehmer test;  fields: exponent,how_far_factored,has_been_pminus1ed
DoubleCheck=11213,0,1     # LL double-check
PRP=1,2,2203,-1,0,0       # PRP;  fields: k,b,n,c[,how_far_factored,tests_saved]
Factor=47000011,50,51     # trial factoring, bit range 50 -> 51
Pminus1=...               # see undoc.txt / readme.txt upstream
```

For TF work units at low bit levels add `UseMaxSieverAllocs=0` (documented in `undoc.txt`), otherwise
the default siever allocation makes small TF ranges take ~700× longer.

**Important for machines that are not 64-bit x86:** the executable reports itself the same way the
upstream Linux64 build does (`Linux64,...`, `PORT=8` / `Wh8` in results) so PrimeNet treats it as a
normal Linux 64-bit client; the CPU-architecture line says `ARM64/AArch64`.  `-DX86_64` is retained
from the upstream makefile only because in the shared sources it means "64-bit build".

---

## Validation (what was actually executed)

### Real ARM64 hardware

Tested on the user's device: **Rockchip RK3528** TV box (model CD1000), Android 9, kernel
**4.19.232 aarch64** (64-bit kernel with a 32-bit Android/Debian userland), **4× Cortex-A53**
(ARMv8, CPU part `0xd03`) up to 2.016 GHz, 2 GB RAM.  The shipped static aarch64 binary runs
unmodified (the 64-bit kernel executes ELF64 even though `uname -m` reports `armv8l` for 32-bit
processes):

* `./mprime -v` → `Mersenne Prime Test Program: Linux64,Untrusted Prime95,v31.7,build 2`
* LL acceptance re-run on real hardware: **14/14** — `M11`/`M67` "is not prime"
  (Res64 `0x6C8` / `677D24EE8AE3B2C2`) and `M127 … M11213` all "is prime!".
  Log: `dist-arm64/tests/realhw-rk3528-ll-results.txt`.
* `mprime -b` benchmark on that device (CPU reported as 2016 MHz, hwloc detected the 4 cores):
  3072K FFT 4 cores/1 worker = **196.4 ms/iter** (5.09 iter/s), 3584K = 249.9 ms/iter,
  4096K = 287.1 ms/iter.  That is roughly two orders of magnitude slower than a modern x86-64
  core for LL/PRP — expected for a 128-bit-SIMD-only port on an in-order Cortex-A53.
  Log: `dist-arm64/tests/realhw-rk3528-bench-and-tf.log`.

### Trial factoring speed on that device (single core; C port of factor64.asm)

| worktodo entry | work | measured time |
|---|---|---|
| `Factor=47000011,50,51` | 1 bit level at 2^50 | **1 s** |
| `Factor=47000011,50,55` | 5 bit levels (2^50→2^55) | **11 s** |
| `Factor=47000011,57,58` | 1 bit level at 2^57 | **44 s** |
| `Factor=47000011,64,65` | 1 bit level at 2^64 (first 2-limb Montgomery level) | **> 244 s**, not finished — extrapolates to ~16–47 h |

TF cost per bit level grows like 2^bits / exponent, and above 2^64 our portable C uses the
documented 2-limb Montgomery path (10–30× slower than the 1-limb path).  Practical GIMPS TF
assignments are 2^70…2^77 at exponents 1.0–1.4e8, i.e. 10^5–10^6 times more work than the
`2^57` case above — **TF on an ARM CPU is not viable**; GIMPS runs TF on GPUs
(mfaktc/mfakto).  Use LL/PRP/P-1 on ARM64, or skip TF work units.  Full log:
`dist-arm64/tests/realhw-rk3528-tf-timing.log`.

### Work-type suitability on real hardware (RK3528, measured)

| Work type | Test | Result on device |
|---|---|---|
| **PRP** (recommended) | `PRP=1,2,216091,-1,0,0` | **53 s**, 0.237 ms/iter (12K FFT), Gerbicz checks 0 errors, `M216091 is a probable prime!`, proof power 5 generated in ~0.5 s (2073 squarings) → best choice |
| LL | `Test=216091,0,1` | **53 s**, 0.237 ms/iter, `M216091 is prime!` — correct, but legacy vs PRP |
| ECM | `ECM=1,2,1277,-1,50000,1000000,1` | stage 1 OK (1.36 s, 1.2 M transforms) then **stage 2 dies silently** (no result) |
| P-1 | `Pminus1=1,2,67,-1,50000,1000000` | stage 1 OK (0.22 s) then fails: `In write_giant, unexpected len == 0 failure` |
| TF | `Factor=...` | minutes → hours → days per bit level at realistic sizes (see TF table above) |

**ECM and P-1 are currently NOT usable on this port**: their stage-2/state paths fail (the scalar
`polymult` implementation documented above as unvalidated is the prime suspect; a state-save
`write_giant len==0` bug also appears for very small moduli).  LL and PRP are the validated paths,
and **PRP is the recommended work type** (same speed as LL, Gerbicz error checking, cheap proof
generation; validated correct on the device).

Practical scale on the RK3528 (from `mprime -b`): 3072K FFT = 196 ms/iter with 4 cores.  A current
GIMPS wavefront PRP (~136M exponent, ~7.5M FFT) would take years on this box; small exponents
(e.g. 216091 took 53 s, ~1M exponent ≈ 20–30 min) are the realistic range.

### Emulated validation (function-correctness, no real hardware needed)

Harness / environment: QEMU 11.1 `qemu-system-aarch64`, `-M virt -cpu cortex-a72 -smp 4`,
Alpine Linux 6.24 aarch64 kernel 6.18.60 as initramfs with the **shipped static mprime binary**
(TCG emulation is function-correct, not speed-representative).

### gwnum FFT math (`_work/gwnum/gwlltest.c`, same converted C sources both targets)

| target | result |
|---|---|
| x86-64 native (mingw) | 12/12 known Mersenne primes (127…11213) → LL residue **0**; p=2047 → `0x8F4A38695034C7E4` (nonzero) |
| aarch64 under QEMU | identical residues, all pass (p=11213: 1.04 s emulated) |

Also verified with mprime's `gwstartnextfft` path, 4 threads, and extended primes up to 216091.
Round-off max 0.1875 (limit 0.5).

### End-to-end mprime on aarch64 (QEMU, static binary)

**Lucas-Lehmer** (`tests/ACC-LL.log`):

```
M11 is not prime. Res64: 00000000000006C8.        <- composite exponent, correct
M67 is not prime. Res64: 677D24EE8AE3B2C2.        <- composite exponent, correct
M127 is prime! ... M11213 is prime!               <- 12/12 known Mersenne primes
```

**PRP with Gerbicz error check** (`tests/ACC-PRP.log`):

```
Starting Gerbicz error-checking PRP test of M2203 using FFT length 112
Gerbicz error check passed at iteration 2116. / 2197.
M2203 is a probable prime!
Starting Gerbicz error-checking PRP test of M67 using Barrett reduction FFT length 32
M67 is not prime.  RES64: 2E99406CF50FC7F1.
```

**Trial factoring** (the C port of `factor64.asm`, `tests/ACC-TF.log`):

```
M1000000009 has a factor: 350544003154897 (TF:44-58)   <- independently predicted by Python
M47000011 no factor from 2^50 to 2^51, Wh8: 00000000  <- independently verified range
```

**P-1** ran on all LL exponents (`M1279 completed P-1, B1=1000, B2=100000` etc.).

Reproduce the mprime tests with the test rig in `../_test/`:
`_test/init-mprime-final.sh` + `_test/payload-mprime/mprime/{worktodo.txt,prime.txt}` and
`_test/mkinitramfs.py`; serial logs are captured in `../_test/logs/`.

---

## What changed vs. the official v31.07b02 tree

`git log`/`git diff` in this directory shows everything (baseline commit = pristine upstream source).

* `gwnum/**` – replaced with the ARM64-capable gwnum 31.7: `_upstream/gwnum_c1` converted assembly
  (`hg/`, `r4/`, `x/` C kernel trees, `xarch.h` NEON/SSE compatibility layer) plus fixes:
  `gwcommon.h` 64-bit mode macro for aarch64; `polymult.c` arch-vs-64bit macro separation
  (`POLYMULT_X86_64`) so the scalar implementation builds on ARM64; `cpuid.c` real Linux/aarch64
  brand/speed/cores/HT/caches/MIDR detection and x86 CPUID helper stubs; `OutputBothRoutine`
  restored (it was hidden behind `GDEBUG` in the converted sources); build matches upstream `make64`
  composition (`-DUSE_SQLITE3` for gwbench's benchmark DB, no `fftsg.c`, no duplicate `mt19937ar`).
  New: `gwnum/Makefile.arm64`; `polymult.c`/`polymult.h` taken from official v31.07.
* `linux-arm64/**` – new platform directory (ports of `linux64/{prime.c,menu.c,os_routines.c,prime.h}`,
  `Makefile.arm64`, `build-arm64-static.{sh,ps1}`, and `factor64_c.c`, the portable C replacement for
  the x86-64 `factor64.asm`, plus a dev-only `factor64_stub.c` for link tests).
* `commonc.c` – 7 guarded lines so the startup line reports `ARM64/AArch64` instead of
  "Not Intel and not AMD" (x86 builds unchanged).
* `security.{c,h}`, `secure5.c` – empty stubs created the same way the upstream makefile does.

Not committed (by design, in `.gitignore`): compiled `.a` libraries and build outputs, which live in
`../_work/` and `../_deps/`.

## License

This is a port of the Prime95 sources; the upstream license applies (see `license.txt` and
`gwnum/LICENSE`-equivalent text inside `gwnum/README.txt`).  If used to find Mersenne primes, GIMPS is
considered the discoverer and the GIMPS prize rules apply.  The ARM64 gwnum conversion is
`sam8641/gwnum_c1`, same license terms.

---

## 10. 2026 round: PrimeNet status, P-1/ECM fixes, ARM instruction-set tuning

### PrimeNet (source-built clients are excluded by the server)
`primenet.c` now sends the `ss` (int) and `sh` (32-char) parameters that PrimeNet v5
requires on every transaction (open-source/untrusted builds send the null values, as the
API specifies).  With that fixed, the server still answers
`pnErrorResult=9: Untrusted program versions currently excluded by PrimeNet` — i.e. a
source-built mprime cannot currently fetch assignments or upload results automatically.
Only trusted (signed) clients can, and none exists for ARM64.  Manual assignments and
manual result submission through the PrimeNet web pages remain possible.
Evidence: `dist-arm64/tests/realhw-rk3528-primenet*.txt`.

### P-1 stage 2 fixed; two bugs remain
* **Fixed:** P-1 stage 2 used to crash because the scalar polymult layout (one double per
  vector, upstream's 32-bit path) does not match gwnum's ARM64 FFT layout (two doubles).
  New `gwnum/arm64_sse2.h` maps the seven `_mm_*_pd` intrinsics polymult uses to GCC
  vector types (compiled to NEON) and ARM64 builds polymult's SSE2 variant; the
  dispatcher now uses it.  Result: `M1277` P-1 runs stage 2 to completion
  (29828 transforms, 0.132 s).
* **Remaining:** `In write_giant, unexpected len == 0 failure` in the P-1 state-save path
  (seen for small exponents M11/M67).  P-1 still "completes" but can report a wrong
  "no factor" (M11 = 2047 = 23*89 with B1=50000 must find 23/89; it reports NF).
* **Remaining:** ECM stage 2 segfaults (SIGSEGV) right after stage-2 initialization;
  ECM stage 1 is correct (M11: found factor 23).

### ARM instruction-set tuning: compiler flags give nothing
Rebuilt gwnum and mprime with `-O3 -mcpu=cortex_a53`: no measurable gain
(LL M216091: 53 s / 0.237 ms-iter at -O2 vs 54 s / 0.242 at -O3; TF 2^57→58: 44 s both).
The converted kernels are already 128-bit vector C that compiles to NEON, and the
in-order Cortex-A53 is instruction/latency limited rather than optimization limited.
Further speedups would need hand-written NEON kernels (est. 1.5–2x for a 2-way TF
Montgomery; ~1.2–1.5x for gwnum FFT), still 1–2 orders of magnitude behind a GPU.

## ECM stage 2 fixed (2026-11-09)

Root cause of the ECM stage-2 SIGSEGV (addr=0x0 inside `mQ_preinit_array`'s std::map
insert): `Dmultiple_data` in `ecm.cpp` did not initialize its `gwdata` and
`do_not_free_val` members, so map-entry copies carried garbage; the optimizer's
bool assumptions turned that into a null-pointer dereference (a -O0/UBSan build
reported "load of value 200, which is not valid for type 'bool'" at the same place).

Fix: initialize `gwdata = NULL; do_not_free_val = TRUE;` in the constructor.
Validated on RK3528: `M1277` ECM curve completes stage 2 (99691 transforms,
0.305 s); ECM on `M67` finds the known factor 193707721 in stage 2; P-1 `M1277`
still completes (regression pass).
