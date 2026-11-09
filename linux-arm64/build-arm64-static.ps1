# ---------------------------------------------------------------------------
# build-arm64-static.ps1 -- reproducible static Linux/aarch64 (musl) mprime
# build, for a Windows host (PowerShell).  The .sh sibling does the same on a
# Unix/git-bash host.
#
# Produces a self-contained ARM64 Linux executable:
#     D:\data\code\p95\_work\mprime\mprime-aarch64   (static ELF64, EM_AARCH64)
#
# Builds, in order:
#   1. ..\gwnum\Makefile.arm64        -> gwnum.a + polymult.a
#      (gwnum 31.7 converted-to-C port, official v31.07 polymult,
#       built with -DUSE_SQLITE3 like upstream Linux builds)
#   2. .\Makefile.arm64 CROSS=1 STATIC=1 -> mprime static link against
#      _deps\arm64\lib\{libhwloc.a,libcurl.a,libgmp.a}
#
# PREREQUISITES
#   * zig 0.15.x (default: ..\..\_tools\zig-x86_64-windows-0.15.1\zig.exe;
#                  override with -Zig or $env:ZIG)
#   * mingw32-make (or use -Make to point at another GNU make)
#   * _deps\arm64\{include,lib} cross-built aarch64-musl GMP, hwloc and curl
#   * Boost headers: _deps\arm64\include\boost (preferred) or the fallback
#     header copy in _work\mprime\boost_1_76_0 (see Makefile.arm64 BOOST_INC)
#
# USAGE (from anywhere):
#   powershell -ExecutionPolicy Bypass -File linux-arm64\build-arm64-static.ps1
#   powershell ... -File ... -Zig C:\zig\zig.exe -Make mingw32-make.exe
# ---------------------------------------------------------------------------
param(
    [string]$Zig  = '',
    [string]$Make = 'mingw32-make',
    [string]$Deps = '../../_deps/arm64',
    [string]$Exe  = '../../_work/mprime/mprime-aarch64'
)

$ErrorActionPreference = 'Stop'
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
Push-Location $here
try {
    $root = (Resolve-Path ../..).Path

    if (-not $Zig) { $Zig = Join-Path $root '_tools\zig-x86_64-windows-0.15.1\zig.exe' }
    if (-not (Test-Path $Zig)) { throw "zig not found: $Zig (pass -Zig C:\path\to\zig.exe)" }
    $zigUnix = (Resolve-Path $Zig).Path -replace '\\', '/'

    if (-not (Get-Command $Make -ErrorAction SilentlyContinue)) {
        throw "'$Make' not on PATH; pass -Make <path to GNU make, e.g. mingw32-make.exe>"
    }

    foreach ($f in @("$Deps/lib/libgmp.a", "$Deps/lib/libhwloc.a", "$Deps/lib/libcurl.a")) {
        if (-not (Test-Path $f)) {
            throw "missing dependency $f (cross-built by the deps-toolchain work item)"
        }
    }

    New-Item -ItemType Directory -Force (Join-Path $root '_work\mprime') | Out-Null

    Write-Host '=== [1/2] gwnum 31.7 aarch64 static libraries ==='
    Push-Location (Join-Path $here '..\gwnum')
    try {
        & $Make -f Makefile.arm64 "ZIG=$zigUnix"
        if ($LASTEXITCODE -ne 0) { throw "gwnum build failed (exit $LASTEXITCODE)" }
    } finally { Pop-Location }
    Get-Item ..\gwnum\gwnum.a, ..\gwnum\polymult.a |
        Select-Object LastWriteTime, Length, Name | Format-Table -AutoSize | Out-Host

    Write-Host '=== [2/2] mprime aarch64 static link ==='
    # The mprime makefile has no header dependency tracking, so start from a
    # clean set of objects every time (guarantees the objects match the current
    # headers, e.g. a re-installed curl/hwloc header tree).
    Get-ChildItem -Path $here -Filter *.o -File | Remove-Item -Force
    & $Make -f Makefile.arm64 CROSS=1 STATIC=1 "ZIG=$zigUnix" "DEPS=$Deps" "EXE=$Exe"
    if ($LASTEXITCODE -ne 0) { throw "mprime build failed (exit $LASTEXITCODE)" }

    if (-not (Test-Path $Exe)) { throw "link did not produce $Exe" }
    $fi = Get-Item $Exe
    Write-Host '=== result ==='
    Write-Host ("  path : " + $fi.FullName)
    Write-Host ("  size : " + $fi.Length + " bytes")
    Write-Host ("  sha256: " + (Get-FileHash $fi.FullName -Algorithm SHA256).Hash)
    $fs = [System.IO.File]::OpenRead($fi.FullName)
    try {
        $buf = New-Object byte[] 20
        $fs.Read($buf, 0, 20) | Out-Null
    } finally { $fs.Close() }
    $em = $buf[18] + 256 * $buf[19]
    Write-Host ("  ELF  : magic=" + (($buf[0..3] | ForEach-Object { $_.ToString('x2') }) -join ' ') +
                " class=" + $buf[4] + " (2=ELF64) data=" + $buf[5] + " (1=LE) e_machine=" + $em + " (183=EM_AARCH64)")
    if ($em -ne 183) { throw "NOT an aarch64 ELF (e_machine=$em)" }
    Write-Host 'OK: static aarch64 mprime built.'
} finally {
    Pop-Location
}
