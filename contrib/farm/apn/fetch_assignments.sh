#!/system/bin/sh
# Fetch ECM assignments from PrimeNet through AutoPrimeNet and hand them to the
# farm (mprime) via worktodo.add, as single-curve lines spread over all workers.
#
# AutoPrimeNet talks to the v5 API in a separate pool directory, using the
# Prime95/MPrime program identity with the untrusted (dummy) security
# parameters: the server rejects the reimplemented trusted-client signature on
# "ga" requests ("Invalid security signature"), while untrusted clients may get
# assignments in the normal security state.  The results themselves are
# submitted separately by submit_apn.sh (trusted-client report-only mode, proven
# to work), keyed by the assignment IDs that mprime keeps in the worktodo lines.
#
# Exponents stay below 131072 so mprime uses the small 5K/8K FFTs: measured on
# the box, a 5K-FFT curve runs at ~224 stage-1 bits/s per worker vs ~90 bits/s
# for the 12K/16K FFT curves (4 workers), i.e. the same work per hour with less
# memory pressure.  See the work-type experiment notes in README.md.
# Environment overrides: MIN_TASKS (default 16 = ~2 days of work), FETCH_N (2), EXP_LO/HI.
MP=/data/adb/mprime
APN=$MP/apn
POOL=$APN/pool
FARM=$MP/farm1
D=/data/local/debian64
LOG=$APN/fetch.log
MIN_TASKS=${MIN_TASKS:-16}
FETCH_N=${FETCH_N:-2}
EXP_LO=${EXP_LO:-100000}
EXP_HI=${EXP_HI:-131072}
log() { echo "$(date '+%F %T') $*" >> "$LOG"; }

# 1) make sure the arm64 chroot and the /data/adb/mprime bind mount exist
/data/local/linux-env/debian64.sh /bin/true 2>/dev/null
if ! mount | grep -q "$D/data/adb/mprime "; then
    mkdir -p "$D/data/adb/mprime"
    mount -o bind /data/adb/mprime "$D/data/adb/mprime" 2>/dev/null || \
        /bin/busybox mount -o bind /data/adb/mprime "$D/data/adb/mprime" 2>/dev/null
fi
if ! mount | grep -q "$D/data/adb/mprime "; then
    log "ERROR: /data/adb/mprime bind mount into $D failed"
    exit 3
fi

# 2) top up only when the farm is running low
tasks=$(cat "$FARM/worktodo.txt" "$FARM/worktodo.add" 2>/dev/null | grep -cE '^(ECM2?|PRP|Pminus1|Pfactor|Test|DoubleCheck|Factor)=')
[ -z "$tasks" ] && tasks=0
if [ "$tasks" -ge "$MIN_TASKS" ]; then
    log "farm has $tasks task lines (>= $MIN_TASKS), no fetch needed"
    exit 0
fi

# 3) one-shot AutoPrimeNet fetch into the pool directory
mkdir -p "$POOL"
[ -f "$POOL/prime.ini" ] || cp -f "$APN/prime.ini" "$POOL/prime.ini"
: > "$POOL/worktodo.txt"
: > "$POOL/results.json.txt"
log "farm has $tasks task lines (< $MIN_TASKS), fetching $FETCH_N assignment(s), exps $EXP_LO-$EXP_HI"
/data/local/linux-env/debian64.sh /usr/bin/python3 "$APN/autoprimenet.py" \
    -w "$POOL" -u vm001 --prime95 --untrusted -T 5 \
    --min-exp "$EXP_LO" --max-exp "$EXP_HI" -n "$FETCH_N" -t 0 \
    -r results.json.txt --no-version-check >> "$LOG" 2>&1
rc=$?
log "AutoPrimeNet fetch rc=$rc"
[ "$rc" -eq 0 ] || exit "$rc"

# 4) split the fetched curves into single-curve lines and append them to
#    worktodo.add (mprime merges that file into worktodo.txt on its own)
/data/local/linux-env/debian64.sh /usr/bin/python3 "$APN/merge_fetch.py" >> "$LOG" 2>&1
rc=$?
log "merge rc=$rc"
exit "$rc"
