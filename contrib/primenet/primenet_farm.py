#!/usr/bin/env python3
"""PrimeNet auto-farm: assign -> run on ARM device -> submit loop.

Drives the manual PrimeNet workflow end to end (source-built mprime is
untrusted for the automatic v5 protocol, so assignments and result upload
go through the logged-in web pages, see primenet_manual.py).

One "pass" does:
  1. make sure the PrimeNet session is alive (login if needed)
  2. if the device has no mprime running: write the next spooled
     assignment(s) to the device worktodo.txt and start mprime -d
  3. pull device results.json.txt; submit any lines not submitted before
  4. if the spool is running low: fetch new assignments from PrimeNet

Run it periodically (cron / Task Scheduler / while-true loop), e.g.:
    set PN_PASS=...
    python primenet_farm.py --host root@rk3528.a0x7fff.site \
        --remote-dir /data/adb/mprime/farm1 --pref 5 --exp-lo 4000000 --exp-hi 9000000

State (spool, submitted-result hashes, session cookies) lives in the
PN_STATE directory (default: directory containing primenet_manual.py).
The PrimeNet password is read from PN_PASS / --password only.
"""

import argparse
import hashlib
import importlib.util
import os
import shlex
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
SPEC = importlib.util.spec_from_file_location(
    "pnm", os.path.join(HERE, "primenet_manual.py"))
pnm = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pnm)


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), msg), flush=True)


def ssh(host, remote_cmd, timeout=60, input_bytes=None):
    cmd = ["ssh", "-o", "ConnectTimeout=10", "-o", "BatchMode=yes",
           host, remote_cmd]
    return subprocess.run(cmd, capture_output=True, timeout=timeout,
                          input=input_bytes)


def scp_to(host, local, remote, timeout=120):
    subprocess.run(["scp", "-o", "ConnectTimeout=10", "-o", "BatchMode=yes",
                    local, "%s:%s" % (host, remote)],
                   capture_output=True, timeout=timeout)


def ensure_login(user, password):
    opener, jar = pnm.make_opener()
    html = opener.open(pnm.BASE + "/manual_result/", timeout=60).read() \
        .decode("utf-8", "replace")
    if 'name="was_logged_in_as"' in html:
        return opener, jar
    log("session expired, logging in again")
    if not password:
        raise SystemExit("PrimeNet session expired and no PN_PASS/--password given")
    if not pnm.cmd_login(opener, jar, user, password):
        raise SystemExit("PrimeNet login failed")
    return pnm.make_opener()


def device_state(host, remote_dir):
    r = ssh(host, "cd %s 2>/dev/null && { ps -ef | grep -v grep | grep -q 'mprime' && echo RUNNING || echo IDLE; wc -l < results.json.txt 2>/dev/null || echo 0; } || echo NODIR" % shlex.quote(remote_dir))
    out = r.stdout.decode().split()
    if not out or out[0] == "NODIR":
        return None
    return {"running": out[0] == "RUNNING",
            "results_lines": int(out[1]) if len(out) > 1 and out[1].isdigit() else 0}


def spool_paths():
    d = pnm.STATE_DIR
    return (os.path.join(d, "farm_spool.txt"),
            os.path.join(d, "farm_submitted.txt"))


def load_spool():
    spool, _ = spool_paths()
    if not os.path.exists(spool):
        return []
    with open(spool, encoding="utf-8") as f:
        return [ln.strip() for ln in f if ln.strip()]


def save_spool(lines):
    spool, _ = spool_paths()
    with open(spool, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + ("\n" if lines else ""))


def top_up_assignments(opener, jar, pref, lo, hi, num, keep_at_least):
    lines = load_spool()
    if len(lines) >= keep_at_least:
        return lines
    want = max(num, keep_at_least - len(lines))
    fields = {"workers": "1", "num_to_get": str(want), "pref": str(pref),
              "exp_lo": str(lo or ""), "exp_hi": str(hi or "")}
    html = pnm.post(opener, "/manual_assignment/", fields)
    got = pnm.ASSIGN_RE.findall(html)
    if got:
        lines.extend(got)
        save_spool(lines)
        log("fetched %d assignment(s): %s" % (len(got), "; ".join(got)))
    else:
        log("PrimeNet returned no assignments for pref=%s range=%s-%s"
            % (pref, lo, hi))
    return lines


def start_device(host, remote_dir, binary, lines, memory_mb):
    spool, _ = spool_paths()
    remote = "%s:%s/" % (host, remote_dir)
    wt_local = os.path.join(pnm.STATE_DIR, "_farm_worktodo.txt")
    with open(wt_local, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    scp_to(host, wt_local, remote + "worktodo.txt")
    # seed a non-interactive prime.txt if the device dir has none
    ssh(host, "cd %s && test -f prime.txt || echo -e 'UsePrimenet=0\\nStressTester=0\\nNumCPUs=1\\nMemory=%d\\nNumWorkers=1\\nWorkPreference=0\\nCoresPerWorker=1' > prime.txt" % (shlex.quote(remote_dir), memory_mb))
    r = ssh(host, "cd %s && test -x mprime || cp %s mprime; chmod +x mprime; nohup ./mprime -d < /dev/null > run.log 2>&1 & echo started" % (shlex.quote(remote_dir), shlex.quote(binary)), timeout=30)
    log("device start: %s" % r.stdout.decode().strip() or r.stderr.decode().strip())


def pull_and_submit(host, remote_dir, opener, jar):
    _, submitted = spool_paths()
    done = set()
    if os.path.exists(submitted):
        with open(submitted, encoding="utf-8") as f:
            done = {ln.split()[0] for ln in f if ln.strip()}
    r = ssh(host, "cat %s/results.json.txt 2>/dev/null" % shlex.quote(remote_dir))
    text = r.stdout.decode("utf-8", "replace")
    new_lines = []
    for ln in text.splitlines():
        ln = ln.strip()
        if not ln:
            continue
        h = hashlib.md5(ln.encode()).hexdigest()
        if h not in done:
            new_lines.append((h, ln))
    if not new_lines:
        return 0
    up = os.path.join(pnm.STATE_DIR, "_farm_submit.json.txt")
    with open(up, "w", encoding="utf-8") as f:
        f.write("\n".join(ln for _, ln in new_lines) + "\n")
    log("submitting %d new result line(s) to PrimeNet" % len(new_lines))
    pnm.cmd_submit(opener, jar, up)
    with open(submitted, "a", encoding="utf-8") as f:
        for h, ln in new_lines:
            f.write("%s %s\n" % (h, ln[:80]))
    return len(new_lines)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--user", default=os.environ.get("PN_USER", "vm001"))
    ap.add_argument("--password", default=os.environ.get("PN_PASS"))
    ap.add_argument("--host", default="root@rk3528.a0x7fff.site")
    ap.add_argument("--remote-dir", default="/data/adb/mprime/farm1")
    ap.add_argument("--binary", default="/data/adb/mprime/mprime-ecm3",
                    help="device path of the mprime binary to copy into the farm dir")
    ap.add_argument("--pref", type=int, default=5)
    ap.add_argument("--exp-lo", type=int)
    ap.add_argument("--exp-hi", type=int)
    ap.add_argument("--num", type=int, default=1, help="assignments per fetch")
    ap.add_argument("--keep", type=int, default=2, help="keep at least this many spooled")
    ap.add_argument("--memory", type=int, default=128, help="device Memory= MB in prime.txt")
    ap.add_argument("--once", action="store_true", help="single pass, no loop")
    args = ap.parse_args()

    while True:
        try:
            opener, jar = ensure_login(args.user, args.password)
            st = device_state(args.host, args.remote_dir)
            if st is None:
                log("device farm dir missing; creating")
                ssh(args.host, "mkdir -p %s" % shlex.quote(args.remote_dir))
                st = {"running": False, "results_lines": 0}
            log("device: %s, results so far: %d" % ("RUNNING" if st["running"] else "IDLE", st["results_lines"]))
            pull_and_submit(args.host, args.remote_dir, opener, jar)
            lines = top_up_assignments(opener, jar, args.pref, args.exp_lo,
                                       args.exp_hi, args.num, args.keep)
            if not st["running"] and lines:
                batch, rest = lines[:1], lines[1:]
                start_device(args.host, args.remote_dir, args.binary,
                             batch, args.memory)
                save_spool(rest)
            elif not st["running"] and not lines:
                log("nothing to run")
        except Exception as e:
            log("ERROR: %r" % e)
        if args.once:
            break
        time.sleep(900)


if __name__ == "__main__":
    main()
