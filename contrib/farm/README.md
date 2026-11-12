# contrib/farm — 7x24 unattended farm operation (RK3528-class Android box)

Self-contained scripts that run **on the device** (no PC required):

- `boot.sh` — boot entry point (idempotent): performance governors (CPU +
  DDR), bind-mounts the farm dir into the debian chroot, starts ddns-go,
  mprime daemon, and the watchdog. Install as an init service via
  `/system/etc/init/farm.rc`:
  ```
  service gimpsfarm /system/bin/sh /data/adb/mprime/boot.sh
      class main
      user root
      group root
      oneshot
  ```
- `farm_watchdog.sh` — forever loop (10 min cycle): re-assert governors,
  restart mprime/ddns-go if dead, submit `results.json.txt` (only truncates
  after a confirmed "Accepted"), top up `worktodo.txt` when a worker's
  section runs low, rotate `run.log` at 50 MB.
- `submit_once.py` — adds the required `checksum` object to each JSON result
  (see ../primenet/README.md) and posts the batch to `/manual_result/`.
  Password: `/data/adb/mprime/.pn_pass` (chmod 600) or `$PN_PASS`.
- `topup_once.py` — fetches ECM assignments (manual form, pref 5) and expands
  each into single-curve `ECM2` lines; rewrites `worktodo.txt` and lets the
  watchdog restart mprime (save files resume automatically).

## Known limitation: PrimeNet WAF vs scripted clients

As of 2026-11, www.mersenne.org's front end **rotates the session cookie on
every request from certain TLS stacks** (observed with Python 3.11 /
OpenSSL 3.0 from the box; Python 3.14 from a Windows host works fine).
Symptom: login succeeds (Set-Cookie issued) but the next request is treated
as anonymous. The watchdog's submission therefore stays best-effort — keep a
PC-side scheduled submission (../primenet) as the primary path until the
server side relaxes this.
