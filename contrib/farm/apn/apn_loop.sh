#!/system/bin/sh
# Keep the farm's PrimeNet result submission running: call submit_apn.sh every
# 30 minutes. Started from boot.sh and supervised by farm_watchdog.sh.
MP=/data/adb/mprime
LOG=$MP/apn/loop.log
log() { echo "$(date '+%F %T') $*" >> "$LOG"; }
log "apn_loop started (pid $$)"
while true; do
    "$MP/apn/submit_apn.sh"
    sleep 1800
done
