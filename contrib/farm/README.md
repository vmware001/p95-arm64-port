# contrib/farm — 7x24 unattended farm operation (RK3528-class Android box)

Self-contained scripts that run **on the device** (no PC required):

- `boot.sh` — boot entry point (idempotent): performance governors (CPU +
  DDR), bind-mounts the farm dir into the debian chroot, starts ddns-go,
  mprime daemon, and the watchdog.
- `farm_watchdog.sh` — forever loop (10 min cycle): re-assert governors,
  restart ddns-go / the AutoPrimeNet loop / mprime if dead, rotate `run.log` at
  50 MB. A singleton lock (`/data/adb/mprime/.watchdog.lock`) keeps exactly one
  copy active: the box's init has been seen to start several copies of the
  service at boot, the extra copies wait and take over if the active one dies.
- `apn/` — **full PrimeNet automation runs on the box itself** (AutoPrimeNet,
  PrimeNet v5 API), in both directions:
  * `apn_loop.sh` — every 30 min calls `submit_apn.sh` and
    `fetch_assignments.sh`; started by `farm_watchdog.sh` and restarted by it
    if it dies (so it comes back automatically after a reboot);
  * `submit_apn.sh` — snapshots `farm1/results.json.txt` into the AutoPrimeNet
    workdir, submits the new lines inside the arm64 Debian chroot
    (`/data/local/debian64`, Python 3.11 + `python3-requests`) and then runs
    `trim_results.py`, which moves the handled lines to
    `farm1/results.submitted.json.txt` and keeps the rest for the next attempt.
    Uses the trusted-client report-only mode (server accepts the signature for
    result reporting);
  * `fetch_assignments.sh` — when the farm worktodo drops below `MIN_TASKS`
    (default 16 task lines = about two days of work at the observed ~8 ECM curves/day), asks AutoPrimeNet for ECM assignments in the
    normal security state, using the Prime95/MPrime program identity with the
    **untrusted** (dummy) security parameters (`--prime95 --untrusted`; the
    reimplemented trusted-client signature is rejected by PrimeNet on `ga`
    requests, see `../primenet/README.md`), restricted to exponents
    `EXP_LO`-`EXP_HI` (100k-400k by default, i.e. the cheap 16K-FFT ECM work);
  * `merge_fetch.py` — splits each fetched assignment into single-curve worktodo
    lines (50 curves -> 50 lines, so the four workers can share it) and appends
    them to `farm1/worktodo.add`, always filling the least-loaded worker first.
    mprime merges `worktodo.add` into `worktodo.txt` on its next save timer
    (<= 10-30 min) and picks the new lines up without a restart;
  * logs: `/data/adb/mprime/apn/{loop,submit,fetch,autoprimenet}.log`.

  Deployment layout on the box (registration details in `../primenet/README.md`):

  ```
  /data/adb/mprime/apn/autoprimenet.py   patched AutoPrimeNet 2.0.2
                                         (../primenet/patches/autoprimenet-2.0.2-*.patch)
  /data/adb/mprime/apn/prime.ini         registration: ComputerGUID + CPU info
                                         (created by the one-time registration command)
  /data/adb/mprime/apn/results_sent.txt  AutoPrimeNet's ledger of handled result lines
  /data/adb/mprime/apn/pool/             AutoPrimeNet workdir for assignment fetching
  /data/adb/mprime/apn/{submit_apn.sh,fetch_assignments.sh,apn_loop.sh,trim_results.py,merge_fetch.py}
  /data/adb/mprime/farm1/worktodo.txt            mprime work spool (single curve per line)
  /data/adb/mprime/farm1/worktodo.add            new assignments from fetch_assignments.sh
  /data/adb/mprime/farm1/results.json.txt        mprime result spool (trimmed after submission)
  /data/adb/mprime/farm1/results.submitted.json.txt   archive of everything handled
  ```

  The arm64 chroot needs `python3-requests` once:
  `debian64.sh /usr/bin/apt-get install -y python3-requests`. The submit/fetch
  scripts re-assert the `/data/adb/mprime` → `/data/local/debian64/data/adb/mprime`
  bind mount on every run.

### Boot autostart

`/system/etc/init/farm.rc` runs the watchdog as an init service
(`class main`, root), the watchdog starts ddns-go, `apn/apn_loop.sh` and mprime
in its first loop:

```
service gimpsfarm /system/bin/sh /data/adb/mprime/farm_watchdog.sh
    class main
    user root
    group root
```

Verified with two cold reboots (2026-11-17): within ~15 s of boot `init.svc.gimpsfarm`
and `init.svc.dropbear` are `running`, the watchdog has restarted ddns-go,
`apn_loop.sh` and mprime, and the first AutoPrimeNet cycle logs
`no new results` (submit) and `farm has N task lines, no fetch needed` (the threshold is now 16)
(fetch). One reboot also verified the `worktodo.add` merge: 150 fetched curves were
folded into `worktodo.txt` (450 task lines) by mprime itself.

### ECM assignments are multi-curve — keep them open

A PrimeNet ECM assignment asks for many curves (CR, typically 50-72). mprime
reports **one result per curve**, so the report must send `d=0` while other
curves of the same assignment are still queued and `d=1` only for the last
curve. The stock AutoPrimeNet 2.0.2 code always sends `d=1` for ECM results, so
PrimeNet closes the assignment when the first curve is reported and every
remaining curve comes back as `PrimeNet error 43: Invalid assignment key`
(observed as `au: no such assignment key` for the heavy assignments). The local
patch `../primenet/patches/autoprimenet-2.0.2-ecm-keep-open.patch` fixes the
`d` flag (factor results still close the assignment).

Consequence for the box (2026-11-17): the four old assignments had already been
closed that way after their first reported curve (M264169, M279863, M7831391,
M10032203), so their ~289 remaining single-curve lines could never be credited
and were dropped. The queue was rebuilt with
`tools/rebalance_worktodo.py` to the two freshly fetched 16K-FFT ECM
assignments (50 curves each, 25 curves per worker); the third fresh assignment
(M289717) was unreserved so the server can hand it out again. After the
rebalance all four workers run the fresh 16K-FFT assignments
(`M244813 curve 1 stage 1 …`), about 1.5-2.5 h per curve per core.

Queue policy: keep the spool small — `fetch_assignments.sh` only tops up below
`MIN_TASKS` (default 16 task lines) and `FETCH_N=2` then adds ~100 curves
(two assignments), i.e. a few days of work at the observed rate, safely inside
the 14-day assignment ETA.

- `submit_once.py` / `topup_once.py` — kept for reference only (old manual-page flow).
  There is no workstation-side submit script any more: submission runs exclusively
  on the box (`apn/` above), so it keeps working with the PC off.

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
