<#
.SYNOPSIS
    Submit Prime95/mprime farm results to PrimeNet through AutoPrimeNet (v5 API, report-only).

.DESCRIPTION
    The mersenne.org *manual results* page is deprecated (James Heinrich, 2026-11-16):
    manual submission is only meant as a last resort, and regular use should go through
    AutoPrimeNet and the PrimeNet v5 API.

    This script implements the farm submission step:
      1. pull /results.json.txt from the ARM box (ssh),
      2. submit the new lines with AutoPrimeNet's Prime95/MPrime report-only mode
         (requires the local patch contrib/primenet/patches/autoprimenet-2.0.2-prime95-report-only.patch
         and a registered computer GUID - see contrib/primenet/README.md),
      3. archive the handled lines locally and remove them from the box so they are
         never submitted twice.

    The script is idempotent: AutoPrimeNet keeps results_sent.txt, and PrimeNet answers
    "error 47 Work no longer needed" for duplicates, which AutoPrimeNet treats as final.

.PARAMETER Box
    ssh destination of the device that runs mprime (default root@tvbox.a0x7fff.site).

.PARAMETER RemoteDir
    mprime working directory on the device (default /data/adb/mprime/farm1).

.PARAMETER ApnDir
    AutoPrimeNet working directory (prime.ini, results_sent.txt, results.json.txt).
    The registration command from contrib/primenet/README.md uses this directory.

.PARAMETER ApnScript
    Path to the patched autoprimenet.py.

.PARAMETER User
    PrimeNet user ID.

.PARAMETER DryRun
    Do everything except the last step: the box's results file is left untouched.

.EXAMPLE
    .\submit_apn.ps1
    .\submit_apn.ps1 -DryRun
#>
[CmdletBinding()]
param(
    [string]$Box = 'root@tvbox.a0x7fff.site',
    [string]$RemoteDir = '/data/adb/mprime/farm1',
    [string]$ApnDir = (Join-Path $PSScriptRoot '..\..\_work\autoprimenet\apn-vm001'),
    [string]$ApnScript = (Join-Path $PSScriptRoot '..\..\_work\autoprimenet\repo\autoprimenet.py'),
    [string]$User = 'vm001',
    [string]$Python = 'python',
    [switch]$DryRun
)

$ErrorActionPreference = 'Stop'
$remoteResults = ($RemoteDir.TrimEnd('/')) + '/results.json.txt'
$localResults = Join-Path $ApnDir 'results.json.txt'
$sentFile = Join-Path $ApnDir 'results_sent.txt'
$logFile = Join-Path $ApnDir 'last_run.log'

if (-not (Test-Path $ApnScript)) { throw "AutoPrimeNet script not found: $ApnScript" }
if (-not (Test-Path $ApnDir)) { New-Item -ItemType Directory -Force -Path $ApnDir | Out-Null }

function Write-AsciiLines([string]$Path, [string[]]$Lines) {
    $text = if ($Lines.Count) { ($Lines -join "`n") + "`n" } else { '' }
    [IO.File]::WriteAllText($Path, $text, [Text.Encoding]::ASCII)
}

Write-Host "[1/4] pulling $Box`:$remoteResults"
$raw = & ssh -o BatchMode=yes -o ConnectTimeout=20 $Box "cat $remoteResults" 2>&1
if ($LASTEXITCODE -ne 0) { throw "ssh failed (rc=$LASTEXITCODE): $raw" }
$lines = @($raw | ForEach-Object { ($_ -replace "`r", '').TrimStart([char]0xFEFF) } | Where-Object { $_.Trim().Length -gt 0 })
Write-AsciiLines $localResults $lines
Write-Host ("      {0} result line(s)" -f $lines.Count)
if ($lines.Count -eq 0) { Write-Host 'Nothing to submit.'; exit 0 }

Write-Host '[2/4] submitting via AutoPrimeNet (v5 API, Prime95/MPrime report-only)'
$apnArgs = @($ApnScript, '-w', $ApnDir, '-u', $User, '--prime95', '-r', 'results.json.txt', '--report-results')
& $Python @apnArgs *>&1 | Tee-Object -FilePath $logFile | Out-Host
$apnRc = $LASTEXITCODE

$credit = 0.0
foreach ($m in Select-String -Path $logFile -Pattern 'CPU credit is ([\d.]+) GHz-days' -AllMatches) {
    foreach ($c in $m.Matches) { $credit += [double]$c.Groups[1].Value }
}
$rejected = (Select-String -Path $logFile -Pattern 'PrimeNet error (?!47)' -AllMatches | Measure-Object).Count
Write-Host ("      AutoPrimeNet rc={0}, credited {1:N3} GHz-days, non-duplicate errors: {2}" -f $apnRc, $credit, $rejected)
if ($apnRc -ne 0) { throw "AutoPrimeNet failed (rc=$apnRc) - see $logFile" }

Write-Host '[3/4] archiving handled lines'
$sent = @()
if (Test-Path $sentFile) {
    $sent = @(Get-Content $sentFile -Encoding ascii | ForEach-Object { ($_ -replace "`r", '').Trim() } | Where-Object { $_.Length -gt 0 })
}
$sentSet = [System.Collections.Generic.HashSet[string]]::new([StringComparer]::Ordinal)
foreach ($s in $sent) { [void]$sentSet.Add($s) }
$handled = @($lines | Where-Object { $sentSet.Contains($_) })
$remaining = @($lines | Where-Object { -not $sentSet.Contains($_) })
if ($handled.Count) {
    $archive = Join-Path $ApnDir ('submitted-' + (Get-Date -Format 'yyyyMMdd-HHmmss') + '.json.txt')
    Write-AsciiLines $archive $handled
    Write-Host ("      {0} handled line(s) archived to {1}" -f $handled.Count, $archive)
} else {
    Write-Host '      no handled lines detected in results_sent.txt'
}

Write-Host '[4/4] updating the farm results file'
if ($remaining.Count -eq $lines.Count) {
    Write-Host '      farm results unchanged'
} elseif ($DryRun) {
    Write-Host ("      DryRun: would leave {0} line(s) on the box" -f $remaining.Count)
} else {
    Write-AsciiLines $localResults $remaining
    $tmp = Join-Path $env:TEMP ('results_push_' + [guid]::NewGuid().ToString('N') + '.txt')
    Copy-Item $localResults $tmp
    try {
        & scp -q $tmp "$Box`:$remoteResults"
        if ($LASTEXITCODE -ne 0) { throw "scp failed (rc=$LASTEXITCODE)" }
        Write-Host ("      farm results updated: {0} line(s) remain" -f $remaining.Count)
    } finally {
        Remove-Item $tmp -Force -ErrorAction SilentlyContinue
    }
}

if ($rejected -gt 0) {
    Write-Warning "There were PrimeNet errors other than duplicate (47) - check $logFile"
    exit 2
}
Write-Host 'Done.'
