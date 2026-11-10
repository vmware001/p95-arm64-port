# PrimeNet manual-workflow automation for the ARM64 mprime build

Source-built mprime is *untrusted* for the PrimeNet v5 automatic protocol
(`pnErrorResult=9 Untrusted program versions currently excluded by PrimeNet`),
and there is no trusted ARM64 client. These scripts drive the **manual**
PrimeNet web workflow (the same pages a logged-in user uses) instead, which
gives automatic assignment fetching and result upload without impersonating a
signed client.

## primenet_manual.py — login / assign / submit

```sh
export PN_PASS=...            # PrimeNet password; never written to disk
python primenet_manual.py --user vm001 login
python primenet_manual.py assign --pref 5 --exp-lo 4000000 --exp-hi 9000000 --out wt.txt
python primenet_manual.py submit --file results.json.txt
```

* `--pref` uses the web-form work codes (150 first-time PRP, 151 PRP DC,
  152 world-record PRP, 4 P-1, 5 ECM, 2 TF, 8 ECM cofactors, ...).
* Upload **must** be `results.json.txt` (JSON); classic text result lines are
  ignored by the server.
* If the requested type has no work in the range, the server returns an
  `ECM2=` line for the range instead.
* Session cookies live in `.pn_cookies.txt` next to the script (or in
  `$PN_STATE` if set); the password is only read from `PN_PASS`/`--password`.

## primenet_farm.py — assign → run on device → submit loop

One pass: keep the session alive, top up a spool of assignments, and if the
ARM device is idle push the next assignment to its `worktodo.txt` and start
`mprime -d`; upload any new `results.json.txt` lines. Designed for cron /
Task Scheduler, e.g. every 15 minutes:

```sh
export PN_PASS=...
python primenet_farm.py --host root@rk3528.a0x7fff.site \
    --remote-dir /data/adb/mprime/farm1 --pref 5 \
    --exp-lo 4000000 --exp-hi 9000000 --keep 2
```

State (spool, submitted-result hashes) lives in `$PN_STATE` (default: this
directory). Add `--once` for a single pass.
