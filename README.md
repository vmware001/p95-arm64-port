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
* `contrib/primenet/` — **manual PrimeNet workflow** scripts (source-built clients are
  excluded from the automatic v5 protocol: `pnErrorResult=9`; these drive the web pages
  instead — login / assign / submit + an auto-farm loop for the RK3528 box)
* `contrib/edtest/` — independent Edwards-arithmetic checker (links only gwnum, verified
  against Python big-int math) that pinned the 16K-FFT `ed_add` bug

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
