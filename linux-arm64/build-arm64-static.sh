#!/bin/sh
# ---------------------------------------------------------------------------
# build-arm64-static.sh -- reproducible static Linux/aarch64 (musl) mprime build
#
# Produces a self-contained `mprime` executable for ARM64 Linux:
#     _work/mprime/mprime-aarch64   (static ELF64, EM_AARCH64)
#
# It builds two components, in order:
#   1. ../gwnum/Makefile.arm64  -> ../gwnum/gwnum.a + ../gwnum/polymult.a
#      (ports of _upstream/gwnum_c1 gwnum 31.7 "asm converted to C", plus the
#       official v31.07 polymult; built with -DUSE_SQLITE3 as upstream does)
#   2. ./Makefile.arm64 CROSS=1 STATIC=1 -> links all mprime objects with
#      _deps/arm64/lib/{libhwloc.a,libcurl.a,libgmp.a} statically.
#
# PREREQUISITES
#   * zig 0.15.x  (default: ../../_tools/zig-x86_64-windows-0.15.1/zig.exe;
#                  set ZIG=/path/to/zig for another install)
#   * GNU make    (on a Windows host via git-bash use MAKE=mingw32-make.exe)
#   * _deps/arm64/{include,lib} cross-built aarch64-musl GMP, hwloc, curl
#     (built by the deps-toolchain work item; see _deps/build-deps.ps1)
#   * Boost headers: _deps/arm64/include/boost, or the fallback header copy in
#     _work/mprime/boost_1_76_0 (see Makefile.arm64 BOOST_INC); or set
#     BOOST_INC=-I/path/to/boost
#
# USAGE (from anywhere):
#   sh linux-arm64/build-arm64-static.sh
#   MAKE=mingw32-make.exe sh linux-arm64/build-arm64-static.sh   # Windows+git-bash
#   ZIG=/usr/local/bin/zig MAKE=make sh linux-arm64/build-arm64-static.sh
#
# All make variables can be overridden from the environment, e.g.
#   DEPS=../../_deps/arm64 EXE=/tmp/mprime JOBS=8
# ---------------------------------------------------------------------------
set -e

HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
cd "$HERE"

ZIG=${ZIG:-"$HERE/../../_tools/zig-x86_64-windows-0.15.1/zig.exe"}
MAKE=${MAKE:-make}
DEPS=${DEPS:-../../_deps/arm64}
EXE=${EXE:-../../_work/mprime/mprime-aarch64}
TARGET=${TARGET:-aarch64-linux-musl}

# Sanity checks with actionable messages.
if [ ! -x "$ZIG" ] && ! command -v "$ZIG" >/dev/null 2>&1; then
	echo "ERROR: zig not found: $ZIG" >&2
	echo "       set ZIG=/path/to/zig" >&2
	exit 1
fi
if ! command -v "$MAKE" >/dev/null 2>&1; then
	echo "ERROR: '$MAKE' not found; set MAKE=mingw32-make.exe (Windows) or MAKE=make" >&2
	exit 1
fi
for f in "$DEPS/lib/libgmp.a" "$DEPS/lib/libhwloc.a" "$DEPS/lib/libcurl.a"; do
	if [ ! -f "$f" ]; then
		echo "ERROR: missing dependency $f" >&2
		echo "       (cross-built by the deps-toolchain work item)" >&2
		exit 1
	fi
done

mkdir -p ../../_work/mprime

echo "=== [1/2] gwnum 31.7 aarch64 static libraries ==="
( cd ../gwnum && "$MAKE" -f Makefile.arm64 ZIG="$ZIG" TARGET="$TARGET" )
ls -l ../gwnum/gwnum.a ../gwnum/polymult.a

# The mprime makefile has no header dependency tracking, so start from a clean
# set of objects every time (guarantees the objects match the current headers,
# e.g. a re-installed curl/hwloc header tree).
rm -f ./*.o
echo "=== [2/2] mprime aarch64 static link ==="
"$MAKE" -f Makefile.arm64 CROSS=1 STATIC=1 \
	ZIG="$ZIG" TARGET="$TARGET" DEPS="$DEPS" EXE="$EXE"

if [ ! -f "$EXE" ]; then
	echo "ERROR: link did not produce $EXE" >&2
	exit 1
fi

echo "=== result ==="
ls -l "$EXE"
# Best-effort ELF verification (readelf is optional).
if command -v readelf >/dev/null 2>&1; then
	readelf -h "$EXE" | grep -E 'Class:|Data:|Machine:|Type:'
else
	# Fallback: check ELF magic + e_machine byte (183 = EM_AARCH64) with od.
	od -An -tx1 -N20 "$EXE" | tr -s ' ' | cut -d' ' -f1-17
	echo "(expect ELF magic 7f 45 4c 46, class 02, little-endian, e_machine b7 at byte 18)"
fi
echo "OK: $EXE"
