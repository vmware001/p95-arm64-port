@echo off
set ZIG_LOCAL_CACHE_DIR=D:\data\code\p95\_work\asm\zigcache\local
set ZIG_GLOBAL_CACHE_DIR=D:\data\code\p95\_work\asm\zigcache\global
set ZIG=D:\data\code\p95\_tools\zig-x86_64-windows-0.15.1\zig.exe
set MAKE=C:\ProgramData\chocolatey\lib\mingw\tools\install\mingw64\bin\mingw32-make.exe
set SRC=D:\data\code\p95\_work\asm
set F=-O3 -ffp-contract=on -mcpu=cortex_a53
set OBJ=%SRC%/obj/obj-base
set LIB=%SRC%/lib/lib-base
set LOG=%SRC%\build_base.log
echo START > %LOG%
if not exist "%SRC%\obj\obj-base" mkdir "%SRC%\obj\obj-base"
if not exist "%SRC%\lib\lib-base" mkdir "%SRC%\lib\lib-base"
cd /d %SRC%\gwnum
"%MAKE%" -j4 -f Makefile.arm64 CC="%ZIG% cc -target aarch64-linux-musl" CXX="%ZIG% c++ -target aarch64-linux-musl" AR="%ZIG% ar" CFLAGS="%F% -Iqd -DUSE_SQLITE3 -I../sqlite-amalgamation-3420000" CXXFLAGS="%F% -Iqd -DUSE_SQLITE3 -I../sqlite-amalgamation-3420000" OBJDIR=%OBJ% LIBDIR=%LIB% >> %LOG% 2>&1
echo MAKE_EXIT=%ERRORLEVEL% >> %LOG%
if not exist "%LIB%\gwnum.a" echo GWNUM_FAIL >> %LOG%
echo DONE >> %LOG%
