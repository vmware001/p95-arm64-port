# PrimeNet manual-workflow automation for the ARM64 mprime build

Source-built mprime is *untrusted* for the PrimeNet v5 automatic protocol
(`pnErrorResult=9 Untrusted program versions currently excluded by PrimeNet`),
and there is no trusted ARM64 client. These scripts drive the **manual**
PrimeNet web workflow (the same pages a logged-in user uses) instead, which
gives automatic assignment fetching and result upload without impersonating a
signed client.

## 2026-11-16: manual results page deprecated — submit through AutoPrimeNet

James Heinrich has asked us to stop using `/manual_result/` for regular
submission (the page itself now says "This page is now considered legacy and
should only be used as a last resort... see AutoPrimeNet for automatic work
management") and to use **AutoPrimeNet** instead. The account's manual results
access was re-enabled after a short risk-control review; the 18 queued ECM
results were submitted that way once and are now clear.

The supported path for this port is AutoPrimeNet's hidden **Prime95/MPrime
report-only mode**, which uses the v5 API with a valid trusted-client signature
(`ss`/`sh`) and permits results whose `"program"` block says
`Prime95 31.7 build 2` (an otherwise untrusted version). Verified 2026-11-17:
a registered AutoPrimeNet instance submitted one of our ECM results and the
server answered `PrimeNet error 47: Work no longer needed / Duplicate result`
— i.e. the API accepted the signature and the result format, and a *new*
result will be credited normally. `pnErrorResult=9` only affects mprime's own
in-client protocol.

The upstream 2.0.2 release validates the program selection without `--prime95`
("Must select a GIMPS program to get assignments for"), so apply the tiny patch
in `patches/autoprimenet-2.0.2-prime95-report-only.patch` first (4 hunks; it
only re-enables the report-only mode that the README still documents). Reported
upstream as <https://github.com/tdulcet/AutoPrimeNet/issues/48> and fixed on the
upstream `dev` branch by <https://github.com/tdulcet/AutoPrimeNet/pull/49>
(`--prime95` on its own selects the Prime95/MPrime program); once a release
contains it the local patch can be dropped.

### Fetching assignments needs the untrusted mode

The reimplemented trusted-client signature is accepted for result reporting
(`t=ar`), but PrimeNet rejects it on assignment requests (`t=ga`) with
`PrimeNet error 9: Access denied / Invalid security signature`.  Untrusted
clients, however, may fetch work, report progress and submit results during the
normal server security state, so apply the second patch
`patches/autoprimenet-2.0.2-untrusted-fetch.patch` and pass the hidden
`--untrusted` option together with `--prime95`: AutoPrimeNet then keeps the
Prime95/MPrime program identity but sends the dummy security parameters
(`ss=19191919`, `sh=ABCD…`).  Verified against the live server: `ga` returns
`PrimeNet success code … Server assigned ECM work` and the assignments are
written as Prime95-style `ECM2=` lines, which our ported mprime understands.

Registration (once, trusted identity is fine) and the two modes:

```sh
python3 autoprimenet.py -w APN -u USER --prime95 --setup
python3 autoprimenet.py -w APN -u USER --prime95 -r results.json.txt --report-results     # upload results
python3 autoprimenet.py -w POOL -u USER --prime95 --untrusted -T 5 -n 2 -t 0 \
    --min-exp 100000 --max-exp 400000 -r results.json.txt --no-version-check             # fetch work
```

```sh
git clone https://github.com/tdulcet/AutoPrimeNet.git
cd AutoPrimeNet
git apply /path/to/p95-arm64-port/contrib/primenet/patches/autoprimenet-2.0.2-prime95-report-only.patch
python -m pip install requests

# 1) one-time registration as a trusted Prime95 client (CPU info of the farm box)
mkdir -p /path/to/apn-vm001 && : > /path/to/apn-vm001/worktodo.txt
python autoprimenet.py -w /path/to/apn-vm001 -u vm001 --prime95 -T 5 \
    -H arm64-rk3528 --processor-model "Rockchip RK3528 (4x Cortex-A53)" \
    --processor-features "ARMv8,ASIMD,NEON" --processor-frequency 2016 \
    --processor-cores 4 --processor-hyperthreads 1 \
    --processor-l1-size 32 --processor-l2-size 512 \
    --memory 1979 --max-memory 768 --hours 24 --timeout 0

# 2) submit the farm's new results (report-only; idempotent, duplicates -> error 47)
python autoprimenet.py -w /path/to/apn-vm001 -u vm001 --prime95 \
    -r results.json.txt --report-results
```

On the box the same thing runs unattended: `contrib/farm/apn/apn_loop.sh` calls
`submit_apn.sh` every 30 minutes inside the arm64 Debian chroot
(`/data/local/debian64`, needs `python3-requests` once), so no workstation is
required. `boot.sh` starts the loop and `farm_watchdog.sh` restarts it if it dies;
log: `/data/adb/mprime/apn/submit.log`.

Assignment fetching for mprime still uses the manual assignment page
(`primenet_manual.py assign` or the web form); only *results* go through the
API. AutoPrimeNet keeps `results_sent.txt`, so lines are never submitted twice.

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
`mprime -d`; upload any new `results.json.txt` lines. Designed for cron / a
shell loop, e.g. every 15 minutes:

```sh
export PN_PASS=...
python primenet_farm.py --host root@YOURBOX.example \
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
