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

ddns-go must be started through `/data/local/linux-env/ddns-go-ctl.sh start`
(the watchdog and `boot.sh` do that): it passes `-c config.yaml -f 300 -dns
223.5.5.5` and appends to `ddns-go.log`. A bare `./ddns-go` starts without the
config file, keeps running but never updates a record — on 2026-11-18 the ISP
re-delegated the IPv6 prefix while ddns-go had been restarted that way, so the
AAAA record stayed on the old prefix and the box became unreachable by name
until ddns-go was restarted through the ctl script.

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

The keep-open logic normally looks for remaining curves in AutoPrimeNet's own
worktodo file, which is empty on the farm: `mprime` owns `farm1/worktodo.txt`
and AutoPrimeNet only reports results. `submit_apn.sh` therefore exports
`APN_ECM_REMAINING_WORKTODO="$MP/farm1/worktodo.txt"` so the remaining curves are
counted from the file `mprime` is actually working from
(`../primenet/patches/autoprimenet-2.0.2-macos-and-farm.patch`). Verified on
2026-11-19 with a real result: 49 curves still queued → the report URL carried
`d=0`, PrimeNet accepted the result (`CPU credit is 0.4858 GHz-days`) and kept
the assignment open; reports with an empty queue (the four test curves before
the fix) carried `d=1`.

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
the 14-day assignment ETA. Exponents are fetched in the 200000-400000 range
(12K/16K FFT) because of the stage-2 memory limit described below.

PrimeNet requests are deliberately rare: `FETCH_MIN_INTERVAL` (default 86400 s)
gates `fetch_assignments.sh` to at most one `ga` attempt per day, even though
`apn_loop.sh` calls it every 30 minutes. On 2026-11-18 a day of frequent
fetch/unreserve/API experiments was followed by PrimeNet's anti-abuse rules
dropping TCP to `209.147.127.187` for hours (ICMP still answered), so the box
now backs off to one attempt a day; a successful attempt queues ~100 curves.

### Work-type experiments (2026-11-17/18)

Goal: find work that yields more credit per hour on the 4x Cortex-A53 box.
GIMPS credit is roughly proportional to the work done, so no work type can
change the FLOP-limited rate much; the measurements confirmed that:

| experiment | result |
|---|---|
| small-exponent ECM (`--min-exp 100000 --max-exp 131072`) | assigns curves with B1≈5.5-5.7M on **5K/8K FFTs**, ~224 stage-1 bits/s/worker vs ~90 bits/s/worker for the 12K/16K work, i.e. about the same credit/hour. **Unusable on this box:** stage 2 allocates ~1.5 GB per worker (`Using 1496MB of memory. D: 51870` in the run log), so mprime answers `Looking for work that uses less memory`, abandons the curve after a whole stage 1 and goes idle. On 2026-11-18 all four workers finished stage 1 of such curves at 02:45-02:56 and then sat idle for ~8 h; those assignments were unreserved. |
| 200k-400k ECM (`--min-exp 200000 --max-exp 400000`) | 12K/16K FFT, B1≈1-2M; stage 2 needs only ~60-560 MB (measured over many curves), so the four workers can serialise stage 2 inside `Memory=768` and finish curves. This is the queue range now. |
| PRP on Mersenne cofactors (work preference 160) | the server accepted the `t=po` preference but answered `t=ga` with ordinary ECM assignments, i.e. it does not hand out PRP-cofactor work to this (untrusted Prime95 v30.19 client) identity; not usable for automation. |
| one worker with `CoresPerWorker=4` (multithreaded single curve, the MT fix path) | a single curve reached ~400 stage-1 bits/s while four independent workers aggregate ~900 bits/s, so the current `NumWorkers=4, CoresPerWorker=1` layout is the faster one. |

Conclusion: keep four single-thread workers and the 12K/16K FFT exponent range;
one small-exponent assignment cannot be completed on a 2 GB box at all.

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
