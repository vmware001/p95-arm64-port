@echo off
set ZIG_LOCAL_CACHE_DIR=D:\data\code\p95\_work\asm\zigcache\local
set ZIG_GLOBAL_CACHE_DIR=D:\data\code\p95\_work\asm\zigcache\global
set ZIG=D:\data\code\p95\_tools\zig-x86_64-windows-0.15.1\zig.exe
set SRC=D:\data\code\p95\_work\asm
set LIB=%SRC%/lib/lib-base
set O=%SRC%/obj/obj-bench
set F=-O3 -ffp-contract=on -mcpu=cortex_a53
set LOG=%SRC%\build_bench.log
echo START > %LOG%
if not exist "%O%" mkdir "%O%"
if not exist "%SRC%\bin" mkdir "%SRC%\bin"
"%ZIG%" cc -target aarch64-linux-musl %F% -I%SRC%\gwnum -c %SRC%\fftbench_asm.c -o %O%\fftbench_asm.o >> %LOG% 2>&1
echo CC1=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl %F% -c %SRC%\prof_asm.c -o %O%\prof_asm.o >> %LOG% 2>&1
echo CC2=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -s -o %SRC%\bin\bench_base %O%\fftbench_asm.o %O%\prof_asm.o %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
echo LINK=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -o %SRC%\bin\bench_nostrip %O%\fftbench_asm.o %O%\prof_asm.o %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
echo LINK_NS=%ERRORLEVEL% >> %LOG%
"%ZIG%" strip -o %SRC%\bin\bench_pad %SRC%\bin\bench_nostrip >> %LOG% 2>&1
echo STRIP=%ERRORLEVEL% >> %LOG%
del %SRC%\bin\bench_pad >> %LOG% 2>&1
dir /b %SRC%\bin >> %LOG% 2>&1
echo DONE >> %LOG%
