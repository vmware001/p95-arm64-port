@echo off
rem Build the multi-variant pass2 experiment harness (one binary, N variants).
set ZIG_LOCAL_CACHE_DIR=D:\data\code\p95\_work\pass2\zigcache\local
set ZIG_GLOBAL_CACHE_DIR=D:\data\code\p95\_work\pass2\zigcache\global
set ZIG=D:\data\code\p95\_tools\zig-x86_64-windows-0.15.1\zig.exe
set SRC=D:\data\code\p95\_work\pass2
set LIB=%SRC%/lib/lib-base
set O=%SRC%/obj/obj-var
set LOG=%SRC%\logs\build_var.log
set F=-O3 -ffp-contract=on -mcpu=cortex_a53
set INC=-I%SRC%/gwnum -I%SRC%/gwnum/qd
echo START > %LOG%
if not exist "%O%" mkdir "%O%"
cd /d %SRC%\gwnum

rem ---- generate restructured copies of r4pass2_8.c ----
cd /d %SRC%
python gen_var8.py loopc gwnum/r4/r4pass2_8_u_loopc.c >> %LOG% 2>&1
echo GEN_LOOPC=%ERRORLEVEL% >> %LOG%
python gen_var8.py all   gwnum/r4/r4pass2_8_u_all.c   >> %LOG% 2>&1
echo GEN_ALL=%ERRORLEVEL% >> %LOG%
python gen_var8.py big   gwnum/r4/r4pass2_8_u_big.c   >> %LOG% 2>&1
echo GEN_BIG=%ERRORLEVEL% >> %LOG%
cd /d %SRC%\gwnum

rem ---- dispatcher ----
"%ZIG%" cc -target aarch64-linux-musl %F% %INC% -c r4/r4vdisp.c -o %O%/r4vdisp.o >> %LOG% 2>&1
echo DISP=%ERRORLEVEL% >> %LOG%

rem ---- 8-level variants ----
"%ZIG%" cc -target aarch64-linux-musl %F% %INC% -DVARIANT=v0 -c r4/r4vwrap_8.c -o %O%/r4v8_v0.o >> %LOG% 2>&1
echo V0=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl %F% %INC% -DVARIANT=v1 -DSRC8="\"r4pass2_8_u_loopc.c\"" -c r4/r4vwrap_8.c -o %O%/r4v8_v1.o >> %LOG% 2>&1
echo V1=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl %F% %INC% -DVARIANT=v2 -DSRC8="\"r4pass2_8_u_all.c\"" -c r4/r4vwrap_8.c -o %O%/r4v8_v2.o >> %LOG% 2>&1
echo V2=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl %F% %INC% -DVARIANT=v3 -DSRC8="\"r4pass2_8_u_big.c\"" -c r4/r4vwrap_8.c -o %O%/r4v8_v3.o >> %LOG% 2>&1
echo V3=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl %F% %INC% -DVARIANT=v4 -DGWNUM_ARM_PREFETCH -DGWNUM_ARM_PREFETCH_LEAD=256 -c r4/r4vwrap_8.c -o %O%/r4v8_v4.o >> %LOG% 2>&1
echo V4=%ERRORLEVEL% >> %LOG%

rem ---- 1280-level variants ----
"%ZIG%" cc -target aarch64-linux-musl %F% %INC% -DVARIANT=v0 -c r4/r4vwrap_1280.c -o %O%/r4v1280_v0.o >> %LOG% 2>&1
echo W0=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl %F% %INC% -DVARIANT=v1 -DGWNUM_ARM_PREFETCH -DGWNUM_ARM_PREFETCH_LEAD=256 -c r4/r4vwrap_1280.c -o %O%/r4v1280_v1.o >> %LOG% 2>&1
echo W1=%ERRORLEVEL% >> %LOG%

rem ---- link harness (variant objects override the archive members) ----
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -s -o %SRC%\bin\bench_var ^
  %SRC%\obj\obj-base\fftbench_asm.o %SRC%\obj\obj-base\prof_asm.o ^
  %O%\r4vdisp.o %O%\r4v8_v0.o %O%\r4v8_v1.o %O%\r4v8_v2.o %O%\r4v8_v3.o %O%\r4v8_v4.o ^
  %O%\r4v1280_v0.o %O%\r4v1280_v1.o ^
  %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
echo LINK=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -o %SRC%\bin\bench_var_nostrip ^
  %SRC%\obj\obj-base\fftbench_asm.o %SRC%\obj\obj-base\prof_asm.o ^
  %O%\r4vdisp.o %O%\r4v8_v0.o %O%\r4v8_v1.o %O%\r4v8_v2.o %O%\r4v8_v3.o %O%\r4v8_v4.o ^
  %O%\r4v1280_v0.o %O%\r4v1280_v1.o ^
  %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
echo LINK_NS=%ERRORLEVEL% >> %LOG%
dir /b %SRC%\bin >> %LOG% 2>&1
echo DONE >> %LOG%
