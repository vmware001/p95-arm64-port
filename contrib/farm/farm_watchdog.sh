#!/system/bin/sh
# GIMPS farm watchdog — keep governors/mprime/ddns-go alive, submit ECM
# results to PrimeNet (manual workflow), top up worktodo when low.
# Runs forever; started by boot.sh. Log: /data/adb/mprime/watchdog.log
PATH=/sbin:/vendor/bin:/system/sbin:/system/bin:/system/xbin
MP=/data/adb/mprime
FARM=$MP/farm1
PN=$MP/pn
PY=/data/local/linux-env/debian.sh
PYARGS="/usr/bin/python3.11"
LOG=$MP/watchdog.log

log() { echo "$(date '+%F %T') $*" >> $LOG; }

while true; do
    # --- governors ---
    for c in 0 1 2 3; do
        g=$(cat /sys/devices/system/cpu/cpu$c/cpufreq/scaling_governor 2>/dev/null)
        [ "$g" != "performance" ] && echo performance > /sys/devices/system/cpu/cpu$c/cpufreq/scaling_governor 2>/dev/null
    done
    g=$(cat /sys/class/devfreq/dmc/governor 2>/dev/null)
    [ "$g" != "performance" ] && echo performance > /sys/class/devfreq/dmc/governor 2>/dev/null

    # --- mprime ---
    if ! pgrep -x mprime >/dev/null 2>&1; then
        log "mprime down; restarting"
        (cd $FARM && nohup ./mprime -d >> run.log 2>&1 &)
        sleep 60
    fi

    # --- ddns-go ---
    if ! pgrep -f '/data/local/ddns-go/ddns-go' >/dev/null 2>&1; then
        log "ddns-go down; restarting"
        (cd /data/local/ddns-go && nohup ./ddns-go >/dev/null 2>&1 &)
    fi

    # --- submit new results ---
    if [ -s $FARM/results.json.txt ]; then
        if $PY $PYARGS $PN/submit_once.py >> $LOG 2>&1; then
            log "submitted one batch"
            : > $FARM/results.json.txt
        else
            log "submit failed"
        fi
    fi

    # --- top up worktodo ---
    $PY $PYARGS $PN/topup_once.py >> $LOG 2>&1

    # --- rotate run.log if > 50 MB ---
    sz=$(stat -c %s $FARM/run.log 2>/dev/null)
    if [ "${sz:-0}" -gt 52428800 ]; then
        tail -c 2000000 $FARM/run.log > $FARM/run.log.tmp
        mv $FARM/run.log.tmp $FARM/run.log
        log "rotated run.log"
    fi

    sleep 600
done
