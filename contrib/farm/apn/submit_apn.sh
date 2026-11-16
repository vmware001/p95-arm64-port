#!/system/bin/sh
# Submit new mprime results to PrimeNet through AutoPrimeNet (v5 API, report-only).
# Runs on the box itself: the arm64 Debian chroot (/data/local/debian64) provides
# Python 3.11 + requests; /data/adb/mprime is bind-mounted into the chroot.
# Idempotent: AutoPrimeNet keeps results_sent.txt and PrimeNet answers error 47
# ("Work no longer needed") for duplicates, which AutoPrimeNet treats as final.
MP=/data/adb/mprime
APN=$MP/apn
D=/data/local/debian64
LOG=$APN/submit.log
log() { echo "$(date '+%F %T') $*" >> "$LOG"; }

# 1) make sure the chroot and the farm bind mount exist (debian64.sh sets up proc/sys/dev)
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

RESULTS=$MP/farm1/results.json.txt
[ -s "$RESULTS" ] || { log "no new results"; exit 0; }

# 2) snapshot the results file into the AutoPrimeNet workdir and submit it
: > "$APN/results.json.txt"
cp -f "$RESULTS" "$APN/results.json.txt"
/data/local/linux-env/debian64.sh /usr/bin/python3 "$APN/autoprimenet.py" \
    -w "$APN" -u vm001 --prime95 -r results.json.txt --report-results >> "$LOG" 2>&1
rc=$?
if [ "$rc" -ne 0 ]; then
    log "AutoPrimeNet failed rc=$rc - results left in place"
    exit "$rc"
fi

# 3) archive handled lines (submitted or definitively rejected) and keep the rest
if [ -s "$APN/results_sent.txt" ]; then
    /data/local/linux-env/debian64.sh /usr/bin/python3 "$APN/trim_results.py" >> "$LOG" 2>&1
    log "trim rc=$?"
fi
exit 0
