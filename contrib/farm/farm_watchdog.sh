#!/system/bin/sh
# GIMPS farm watchdog — NETWORK-FREE build.
#
# Keeps governors / mprime / ddns-go alive and rotates the log.  It
# deliberately does NOT talk to PrimeNet and does NOT top up worktodo:
# the vm001 account is under risk-control review, so *no* request to
# mersenne.org may originate from the box until the flag is lifted.
# Results are submitted later, in one manual batch, from the workstation.
#
# Runs forever; started by boot.sh (or as an init service, see farm.rc).
PATH=/sbin:/vendor/bin:/system/sbin:/system/bin:/system/xbin
MP=/data/adb/mprime
FARM=$MP/farm1
LOG=$MP/watchdog.log

log() { echo "$(date '+%F %T') $*" >> $LOG; }

while true; do
    # --- governors (this box resets them occasionally) ---
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

    # --- ddns-go (only so the workstation can still reach this box over ssh) ---
    if ! pgrep -f '/data/local/ddns-go/ddns-go' >/dev/null 2>&1; then
        log "ddns-go down; restarting"
        (cd /data/local/ddns-go && nohup ./ddns-go >/dev/null 2>&1 &)
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
