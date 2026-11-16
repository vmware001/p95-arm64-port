@echo off
set ZIG_LOCAL_CACHE_DIR=D:\data\code\p95\_work\asm\zigcache\local
set ZIG_GLOBAL_CACHE_DIR=D:\data\code\p95\_work\asm\zigcache\global
set ZIG=D:\data\code\p95\_tools\zig-x86_64-windows-0.15.1\zig.exe
set MAKE=C:\ProgramData\chocolatey\lib\mingw\tools\install\mingw64\bin\mingw32-make.exe
set SRC=D:\data\code\p95\_work\asm
set F=-O3 -ffp-contract=on -mcpu=cortex_a53
set OBJ=%SRC%/obj/obj-asm
set LIB=%SRC%/lib/lib-asm
set LOG=%SRC%\build_asm.log
echo START > %LOG%
if not exist "%SRC%\obj\obj-asm" mkdir "%SRC%\obj\obj-asm"
if not exist "%SRC%\lib\lib-asm" mkdir "%SRC%\lib\lib-asm"
cd /d %SRC%\gwnum
"%MAKE%" -j4 -f Makefile.arm64 CC="%ZIG% cc -target aarch64-linux-musl" CXX="%ZIG% c++ -target aarch64-linux-musl" AR="%ZIG% ar" CFLAGS="%F% -DGWNUM_ASM_PASS1_COPY -Iqd -DUSE_SQLITE3 -I../sqlite-amalgamation-3420000" CXXFLAGS="%F% -DGWNUM_ASM_PASS1_COPY -Iqd -DUSE_SQLITE3 -I../sqlite-amalgamation-3420000" OBJDIR=%OBJ% LIBDIR=%LIB% >> %LOG% 2>&1
echo MAKE_EXIT=%ERRORLEVEL% >> %LOG%
if not exist "%LIB%\gwnum.a" echo GWNUM_FAIL >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl %F% -I%SRC%\gwnum -c %SRC%\fftbench_asm.c -o %SRC%\obj\obj-asm\fftbench_asm.o >> %LOG% 2>&1
"%ZIG%" cc -target aarch64-linux-musl %F% -c %SRC%\prof_asm.c -o %SRC%\obj\obj-asm\prof_asm.o >> %LOG% 2>&1
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -s -o %SRC%\bin\bench_asm %SRC%\obj\obj-asm\fftbench_asm.o %SRC%\obj\obj-asm\prof_asm.o %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
echo LINK=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -o %SRC%\bin\bench_asm_nostrip %SRC%\obj\obj-asm\fftbench_asm.o %SRC%\obj\obj-asm\prof_asm.o %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
"%ZIG%" cc -target aarch64-linux-musl -O2 -static -no-pie -s -o %SRC%\bin\copytest %SRC%\copytest.c %OBJ%\r4pass1_asm.o -lm >> %LOG% 2>&1
echo COPYTEST_LINK=%ERRORLEVEL% >> %LOG%
dir /b %SRC%\bin >> %LOG% 2>&1
echo DONE >> %LOG%
