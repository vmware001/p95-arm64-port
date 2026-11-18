#!/bin/sh
# ===========================================================================
# build-macos-arm64.sh -- build native mprime for Apple Silicon (arm64 macOS)
#
# Run this on an Apple Silicon Mac (M1/M2/M3/M4) with Xcode command line tools:
#
#   brew install gmp hwloc boost curl
#   ./macosxarm64/build-macos-arm64.sh
#
# Produces macosxarm64/mprime (Mach-O arm64, no Rosetta).
#
# Notes:
#  * gwnum is the gwnum 31.7 "assembly converted to C" ARM64 port from
#    ../gwnum (same kernels as the Linux arm64 build in this tree).
#  * APPLE_CPU tunes for the core family (default apple-m1); try apple-m2,
#    apple-m3, apple-m4 or APPLE_CPU= for a generic armv8 build:
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
