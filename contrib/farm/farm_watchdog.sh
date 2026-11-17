#!/system/bin/sh
# GIMPS farm watchdog.
#
# Keeps governors / mprime / ddns-go / the AutoPrimeNet submit loop alive and
# rotates the log.  It does NOT submit anything itself: result submission is
# done by apn/apn_loop.sh (AutoPrimeNet, PrimeNet v5 API, report-only), which
# this watchdog restarts if it ever dies.
#
# Runs forever; started by boot.sh (or as an init service, see farm.rc).
PATH=/sbin:/vendor/bin:/system/sbin:/system/bin:/system/xbin
MP=/data/adb/mprime
FARM=$MP/farm1
LOG=$MP/watchdog.log

log() { echo "$(date '+%F %T') $*" >> $LOG; }

# --- singleton guard ---------------------------------------------------------
# The box's init has been observed to start several copies of this service at
# boot.  Only one copy may act; the others wait and take over when the active
# owner dies.  The lock is a directory (mkdir is atomic).  Waiting copies never
# exit, so init never enters a restart loop.
LOCK=$MP/.watchdog.lock
while true; do
    if mkdir "$LOCK" 2>/dev/null; then
        echo $$ > "$LOCK/pid"
        break
    fi
    owner=$(cat "$LOCK/pid" 2>/dev/null)
    if [ -z "$owner" ]; then
        sleep 5
        continue
    fi
    if kill -0 "$owner" 2>/dev/null; then
        sleep 60
    else
        rm -rf "$LOCK"
        sleep 5
    fi
done

while true; do
    # --- governors (this box resets them occasionally) ---
    for c in 0 1 2 3; do
        g=$(cat /sys/devices/system/cpu/cpu$c/cpufreq/scaling_governor 2>/dev/null)
        [ "$g" != "performance" ] && echo performance > /sys/devices/system/cpu/cpu$c/cpufreq/scaling_governor 2>/dev/null
    done
    g=$(cat /sys/class/devfreq/dmc/governor 2>/dev/null)
    [ "$g" != "performance" ] && echo performance > /sys/class/devfreq/dmc/governor 2>/dev/null

    # --- ddns-go first, so the workstation can reach this box again as soon as
    #     possible after a reboot (only ddns-go is needed for remote access) ---
    if ! pgrep -x ddns-go >/dev/null 2>&1; then
        log "ddns-go down; restarting"
        (cd /data/local/ddns-go && nohup ./ddns-go >/dev/null 2>&1 &)
    fi

    # --- AutoPrimeNet loop (report results + fetch assignments via the v5 API) ---
    if [ -x "$MP/apn/apn_loop.sh" ] && ! pgrep -f "$MP/apn/apn_loop.sh" >/dev/null 2>&1; then
        log "apn_loop down; restarting"
        nohup $MP/apn/apn_loop.sh >/dev/null 2>&1 &
    fi

    # --- mprime ---
    if ! pgrep -x mprime >/dev/null 2>&1; then
        log "mprime down; restarting"
        (cd $FARM && nohup ./mprime -d >> run.log 2>&1 &)
    fi

    # --- rotate run.log if > 50 MB ---
    sz=$(stat -c %s $FARM/run.log 2>/dev/null)
    if [ "${sz:-0}" -gt 52428800 ]; then
        tail -c 2000000 $FARM/run.log > $FARM/run.log.tmp
        mv $FARM/run.log.tmp $FARM/run.log
        log "rotated run.log"
    fi

    sleep 600
done
