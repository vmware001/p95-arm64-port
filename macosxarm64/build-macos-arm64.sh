#!/bin/sh
# ===========================================================================
# build-macos-arm64.sh -- build native mprime for Apple Silicon (arm64 macOS)
#
# Run this on an Apple Silicon Mac (M1/M2/M3/M4/M5/M6) with Xcode command line
# tools:
#
#   brew install gmp hwloc boost curl
#   ./macosxarm64/build-macos-arm64.sh
#
# Produces macosxarm64/mprime (Mach-O arm64, no Rosetta).
#
# IMPORTANT -- hwloc must contain the Darwin perflevel fix:
#  * Stock hwloc 2.11.x / 2.15.0 hangs forever in hwloc_topology_load() on
#    Apple Silicon whose perflevel cpusets are not a multiple of the L2 group
#    width (observed on Apple M6, macOS 27.2).  mprime then starts at 100% CPU
#    and never prints anything.
#  * Apply macosxarm64/patches/hwloc-2.15.0-darwin-perflevel-cache-fix.patch
#    and build hwloc from source, then point HWLOC_PREFIX at it:
#        HWLOC_PREFIX=/path/to/patched/hwloc ./macosxarm64/build-macos-arm64.sh
#  * A pure from-source build (no Homebrew at all) is also supported: build
#    gmp/hwloc/boost into one prefix and set GMP_PREFIX/HWLOC_PREFIX/
#    BOOST_PREFIX.  See macosxarm64/README.md.
#
# Notes:
#  * gwnum is the gwnum 31.7 "assembly converted to C" ARM64 port from
#    ../gwnum (same kernels as the Linux arm64 build in this tree).
#  * APPLE_CPU tunes for the core family (default apple-m1); try apple-m2,
#    apple-m3, apple-m4, apple-m5, apple-m6 or APPLE_CPU= for a generic armv8
#    build:
#        APPLE_CPU=apple-m3 ./macosxarm64/build-macos-arm64.sh
#  * Everything (CC/CPP/CFLAGS/GMP_PREFIX/...) can be overridden from the
#    environment; see the makefile.
# ===========================================================================
set -e

HERE=$(cd "$(dirname "$0")" && pwd)
cd "$HERE"

echo "=== host ==="
uname -a
sysctl -n machdep.cpu.brand_string 2>/dev/null || true
echo "clang: $(clang --version | head -1)"

echo "=== dependencies (Homebrew) ==="
for f in gmp hwloc boost curl; do
	if brew list --versions "$f" >/dev/null 2>&1; then
		echo "  $f: $(brew list --versions "$f")"
	else
		echo "  WARNING: $f not installed; run: brew install $f" >&2
	fi
done

echo "=== build ==="
JOBS=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)
make -j"$JOBS" APPLE_CPU="${APPLE_CPU:-apple-m1}"

echo "=== verify ==="
ls -l mprime
file mprime
otool -hv mprime | tail -3
echo "=== version ==="
./mprime -v || true

echo "OK: $HERE/mprime"
