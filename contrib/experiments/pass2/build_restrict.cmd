@echo off
rem Build variant v5: pass2 TU compiled with __restrict butterfly pointers in r4.h
rem (lets the compiler hoist the 8 interleaved loads above the 8 stores inside
rem each butterfly - the arithmetic itself is untouched).
set ZIG_LOCAL_CACHE_DIR=D:\data\code\p95\_work\pass2\zigcache\local
set ZIG_GLOBAL_CACHE_DIR=D:\data\code\p95\_work\pass2\zigcache\global
set ZIG=D:\data\code\p95\_tools\zig-x86_64-windows-0.15.1\zig.exe
set SRC=D:\data\code\p95\_work\pass2
set LIB=%SRC%/lib/lib-base
set O=%SRC%/obj/obj-var
set LOG=%SRC%\logs\build_restrict.log
set F=-O3 -ffp-contract=on -mcpu=cortex_a53
set INC=-I%SRC%/gwnum -I%SRC%/gwnum/qd
echo START > %LOG%
cd /d %SRC%\gwnum
"%ZIG%" cc -target aarch64-linux-musl %F% %INC% -c r4/r4vdisp.c -o %O%/r4vdisp.o >> %LOG% 2>&1
echo DISP=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl %F% %INC% -DVARIANT=v5 -c r4/r4vwrap_8.c -o %O%/r4v8_v5.o >> %LOG% 2>&1
echo V5=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -s -o %SRC%\bin\bench_var ^
  %SRC%\obj\obj-base\fftbench_asm.o %SRC%\obj\obj-base\prof_asm.o ^
  %O%\r4vdisp.o %O%\r4v8_v0.o %O%\r4v8_v1.o %O%\r4v8_v2.o %O%\r4v8_v3.o %O%\r4v8_v4.o %O%\r4v8_v5.o ^
  %O%\r4v1280_v0.o %O%\r4v1280_v1.o ^
  %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
echo LINK=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -o %SRC%\bin\bench_var_nostrip ^
  %SRC%\obj\obj-base\fftbench_asm.o %SRC%\obj\obj-base\prof_asm.o ^
  %O%\r4vdisp.o %O%\r4v8_v0.o %O%\r4v8_v1.o %O%\r4v8_v2.o %O%\r4v8_v3.o %O%\r4v8_v4.o %O%\r4v8_v5.o ^
  %O%\r4v1280_v0.o %O%\r4v1280_v1.o ^
  %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
echo LINK_NS=%ERRORLEVEL% >> %LOG%
echo DONE >> %LOG%
