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

## add_checksum.py — REQUIRED since 2026 for ECM submission

The `/manual_result/` form rejects ECM JSON without a top-level `checksum`
object: *"Error: Missing checksum"*. Source builds emit
`security-code: 00000000` (the proprietary `security.h` is not shipped), and
that field is not accepted as a checksum.

```sh
python add_checksum.py results.json.txt out.json
python primenet_manual.py --user vm001 submit --file out.json
```

It appends `"checksum":{"version":1,"checksum":"<CRC32-8hex>"}` computed over
`{exp};ECM;;{curves};{b1};{b2};;;;{prog};{ver};;;ARM64 Linux;{timestamp}`
(the open mfaktc TF-checksum layout adapted for ECM). Empirically
(2026-11-12) the server only checks the field's *presence* for ECM — the CRC
value is not re-validated — but it is computed correctly anyway.

## Hard-won lessons (RK3528 farm bring-up, 2026-11)

- ECM results are **batched per worktodo line**: `ECM2=...,N` produces ONE
  result after all N curves. Split into single-curve lines for per-curve
  reporting.
- The server **deletes the assignment** on the first reported curve of an AID;
  further same-AID curves still credit GHz-days but without assignment
  bookkeeping. Keep locally tracking what remains to run.
- **B1 must equal the assignment's bound** or the result is dropped
  (silently). B2 is free — report the memory-adaptive actual; larger earns
  more credit.
- Fetch result files **byte-exact**: a UTF-16LE redirect garbles the JSON and
  the server then silently ignores the whole batch with an empty
  "Done processing:". Verify the first bytes are `{"status"`.
- On duplicate submission the server may double-credit (no robust dedupe for
  NF-ECM) — track submitted hashes in `$PN_STATE` (`submitted.json`).
