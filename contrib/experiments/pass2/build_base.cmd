@echo off
rem Baseline build (identical flags to production mc4) in private dir.
set ZIG_LOCAL_CACHE_DIR=D:\data\code\p95\_work\pass2\zigcache\local
set ZIG_GLOBAL_CACHE_DIR=D:\data\code\p95\_work\pass2\zigcache\global
set ZIG=D:\data\code\p95\_tools\zig-x86_64-windows-0.15.1\zig.exe
set MAKE=C:\ProgramData\chocolatey\lib\mingw\tools\install\mingw64\bin\mingw32-make.exe
set SRC=D:\data\code\p95\_work\pass2
set F=-O3 -ffp-contract=on -mcpu=cortex_a53
set OBJ=%SRC%/obj/obj-base
set LIB=%SRC%/lib/lib-base
set LOG=%SRC%\logs\build_base.log
echo START > %LOG%
if not exist "%SRC%\obj\obj-base" mkdir "%SRC%\obj\obj-base"
if not exist "%SRC%\lib\lib-base" mkdir "%SRC%\lib\lib-base"
if not exist "%SRC%\bin" mkdir "%SRC%\bin"
cd /d %SRC%\gwnum
"%MAKE%" -j4 -f Makefile.arm64 CC="%ZIG% cc -target aarch64-linux-musl" CXX="%ZIG% c++ -target aarch64-linux-musl" AR="%ZIG% ar" CFLAGS="%F% -Iqd -DUSE_SQLITE3 -I../sqlite-amalgamation-3420000" CXXFLAGS="%F% -Iqd -DUSE_SQLITE3 -I../sqlite-amalgamation-3420000" OBJDIR=%OBJ% LIBDIR=%LIB% >> %LOG% 2>&1
echo MAKE_EXIT=%ERRORLEVEL% >> %LOG%
if not exist "%LIB%\gwnum.a" echo GWNUM_FAIL >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl %F% -Iqd -DUSE_SQLITE3 -I../sqlite-amalgamation-3420000 -c polymult.c -o %OBJ%/polymult.o >> %LOG% 2>&1
"%ZIG%" cc -target aarch64-linux-musl %F% -Iqd -DSSE2 -c polymult.c -o %OBJ%/polymult_sse2.o >> %LOG% 2>&1
"%ZIG%" ar -rs %LIB%/polymult.a %OBJ%/polymult.o %OBJ%/polymult_sse2.o >> %LOG% 2>&1
echo POLYMULT=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl %F% -I%SRC%\gwnum -c %SRC%\fftbench_asm.c -o %OBJ%\fftbench_asm.o >> %LOG% 2>&1
"%ZIG%" cc -target aarch64-linux-musl %F% -c %SRC%\prof_asm.c -o %OBJ%\prof_asm.o >> %LOG% 2>&1
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -s -o %SRC%\bin\bench_base %SRC%\obj\obj-base\fftbench_asm.o %SRC%\obj\obj-base\prof_asm.o %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
echo LINK=%ERRORLEVEL% >> %LOG%
"%ZIG%" cc -target aarch64-linux-musl -static -no-pie -o %SRC%\bin\bench_base_nostrip %SRC%\obj\obj-base\fftbench_asm.o %SRC%\obj\obj-base\prof_asm.o %LIB%/gwnum.a %LIB%/polymult.a -lm -lpthread >> %LOG% 2>&1
echo LINK_NS=%ERRORLEVEL% >> %LOG%
dir /b %SRC%\bin >> %LOG% 2>&1
echo DONE >> %LOG%
