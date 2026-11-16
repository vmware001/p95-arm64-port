@echo off
rem Build profiling binary + microbenchmark (separate from timing binaries).
set ZIG_LOCAL_CACHE_DIR=D:\data\code\p95\_work\pass2\zigcache\local
set ZIG_GLOBAL_CACHE_DIR=D:\data\code\p95\_work\pass2\zigcache\global
set ZIG=D:\data\code\p95\_tools\zig-x86_64-windows-0.15.1\zig.exe
set SRC=D:\data\code\p95\_work\pass2
set LIB=%SRC%/lib/lib-base
set O=%SRC%/obj/obj-misc
set F=-O3 -ffp-contract=on -mcpu=cortex_a53
set LOG=%SRC%\logs\build_misc.log
echo START > %LOG%
if not exist "%O%" mkdir "%O%"
"%ZIG%" cc -target aarch64-linux-musl -O2 -g -mcpu=cortex_a53 -c %SRC%\prof2.c -o %O%\prof2.o >> %LOG% 2>&1
echo CC_PROF=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -O3 -DNDEBUG -mcpu=cortex_a53 -c %SRC%\ubench.c -o %O%\ubench.o >> %LOG% 2>&1
echo CC_UB=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -s -o %SRC%\bin\ubench %O%\ubench.o >> %LOG% 2>&1
echo LINK_UB=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -O3 -DNDEBUG -mcpu=cortex_a53 -c %SRC%\ubench2.c -o %O%\ubench2.o >> %LOG% 2>&1
echo CC_UB2=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -s -o %SRC%\bin\ubench2 %O%\ubench2.o >> %LOG% 2>&1
echo LINK_UB2=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -O2 -g -mcpu=cortex_a53 -I%SRC%\gwnum -c %SRC%\fftbench_asm.c -o %O%\fftbench_asm_g.o >> %LOG% 2>&1
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -g -o %SRC%\bin\bench_prof2 %O%\fftbench_asm_g.o %O%\prof2.o %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
echo LINK_PROF=%ERRORLEVEL% >> %LOG%
dir /b %SRC%\bin >> %LOG% 2>&1
echo DONE >> %LOG%
