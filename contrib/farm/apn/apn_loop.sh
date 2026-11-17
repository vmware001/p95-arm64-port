#!/system/bin/sh
# Keep the farm's PrimeNet automation running, every 30 minutes:
#   1. submit any new results through AutoPrimeNet (report-only, trusted client),
#   2. fetch new ECM assignments through AutoPrimeNet when the farm runs low.
# Started by boot.sh / farm_watchdog.sh (it is restarted automatically if it
# ever dies).
MP=/data/adb/mprime
LOG=$MP/apn/loop.log
log() { echo "$(date '+%F %T') $*" >> "$LOG"; }
log "apn_loop started (pid $$)"
while true; do
    "$MP/apn/submit_apn.sh"
    "$MP/apn/fetch_assignments.sh"
    sleep 1800
done
