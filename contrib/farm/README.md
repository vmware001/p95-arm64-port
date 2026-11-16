# contrib/farm — 7x24 unattended farm operation (RK3528-class Android box)

Self-contained scripts that run **on the device** (no PC required):

- `boot.sh` — boot entry point (idempotent): performance governors (CPU +
  DDR), bind-mounts the farm dir into the debian chroot, starts ddns-go,
  mprime daemon, and the watchdog. Install as an init service via
  `/system/etc/init/farm.rc`:
  ```
  service gimpsfarm /system/bin/sh /data/adb/mprime/farm_watchdog.sh
      class main
      user root
      group root
  ```
- `farm_watchdog.sh` — forever loop (10 min cycle): re-assert governors,
  restart mprime/ddns-go if dead, rotate `run.log` at 50 MB. **Network-free**
  (see "account lockdown" below).
- `apn/` — **PrimeNet result submission runs on the box itself** (AutoPrimeNet,
  PrimeNet v5 API, report-only):
  * `apn_loop.sh` — every 30 min calls `submit_apn.sh`; started by `boot.sh`
    (step 6) and restarted by `farm_watchdog.sh` if it dies;
  * `submit_apn.sh` — snapshots `farm1/results.json.txt` into the AutoPrimeNet
    workdir, submits the new lines inside the arm64 Debian chroot
    (`/data/local/debian64`, Python 3.11 + `python3-requests`) and then runs
    `trim_results.py`, which moves the handled lines to
    `farm1/results.submitted.json.txt` and keeps the rest for the next attempt;
  * logs: `/data/adb/mprime/apn/{loop,submit,autoprimenet}.log`.

  Deployment layout on the box (registration details in `../primenet/README.md`):

  ```
  /data/adb/mprime/apn/autoprimenet.py   patched AutoPrimeNet 2.0.2
                                         (../primenet/patches/autoprimenet-2.0.2-prime95-report-only.patch)
  /data/adb/mprime/apn/prime.ini         registration: ComputerGUID + CPU info
                                         (created by the one-time registration command)
  /data/adb/mprime/apn/results_sent.txt  AutoPrimeNet's ledger of handled lines
  /data/adb/mprime/apn/{submit_apn.sh,apn_loop.sh,trim_results.py}
  /data/adb/mprime/farm1/results.json.txt             mprime result spool (trimmed after submission)
  /data/adb/mprime/farm1/results.submitted.json.txt   archive of everything handled
  ```

  The arm64 chroot needs `python3-requests` once:
  `debian64.sh /usr/bin/apt-get install -y python3-requests`. `boot.sh` sets up the
  `/data/adb/mprime` → `/data/local/debian64/data/adb/mprime` bind mount; the submit
  script re-asserts it on every run.
- `submit_once.py` / `topup_once.py` — kept for reference only (old manual-page flow).
- `submit_apn.ps1` — **workstation fallback** (same AutoPrimeNet submission from a PC;
  pulls/pushes the results file over ssh). The box-side `apn/` loop is the primary
  path and needs no PC.

## Run 4 workers, not 3: Memory must be big enough

Symptom (mprime 31.7, 4x A53 / 2 GB box, `Memory=576`, NumWorkers=4):

```
[Worker #4] Looking for work that uses less memory.     (x7177)
[Worker #4] No work to do at the present time.  Waiting.
[Worker #4] Exceeded limit on number of workers that can use lots of memory.
```

and `/proc/stat` shows only 3 of 4 cores busy — Worker #4 parks on a futex and
*never* retries, because once `MEM_RESTART_FLAGS[thread]` is set the worker
skips every worktodo line (`is_worker_waiting_for_more_memory()` in
`commonb.c`). With `Memory=576` the per-worker allowance is 144 MB, which the
ECM stage-2 planner rejects.

Fix: raise the global memory budget so the per-worker share covers stage 2,
e.g. `Memory=768` on a 2 GB box (per-worker allowance becomes ~202 MB, all
four workers start). Measured effect on the box: 73 % -> 97 % total CPU
utilisation (+33 % throughput). Monitor `MemAvailable` afterwards; zram
provides a small swap cushion.

A restart with `kill -TERM` is safe: mprime checkpoints the current
iteration on exit, so no measurable work is lost.

## Account lockdown (2026-11)

The vm001 account was flagged by PrimeNet ("suspicious activity") after a
burst of scripted manual-submission requests. While the review is pending the
box must not issue *any* mersenne.org request. Note that `chmod 000` does not
stop a root-run watchdog — the scripts were still being executed every
10 minutes. The lockdown therefore is:

1. `farm_watchdog.sh` replaced by the network-free version above;
2. `submit_once.py` / `topup_once.py` renamed to `*.py.disabled`;
3. `results.json.txt` left untouched (never truncated without an "Accepted").

## Known limitation: PrimeNet WAF vs scripted clients

www.mersenne.org's front end **rotates the session cookie on every request
from certain TLS stacks** (observed with Python 3.11 / OpenSSL 3.0 from the
box; Python 3.14 from a Windows host works). Symptom: login succeeds
(Set-Cookie issued) but the next request is treated as anonymous. This is why
manual browser submission from a workstation is the primary path.
