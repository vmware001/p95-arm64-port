#!/system/bin/sh
# GIMPS ECM farm boot script — runs as root at boot via
# /system/etc/init/farm.rc (or manually). Idempotent.
PATH=/sbin:/vendor/bin:/system/sbin:/system/bin:/system/xbin
MP=/data/adb/mprime
FARM=$MP/farm1

# 1) CPU + DDR governors (this box resets them occasionally; watchdog re-asserts too)
for c in 0 1 2 3; do
    echo performance > /sys/devices/system/cpu/cpu$c/cpufreq/scaling_governor 2>/dev/null
done
echo performance > /sys/class/devfreq/dmc/governor 2>/dev/null

# 2) bind the farm dir into the armhf debian chroot (linux-env wrappers chroot,
#    so /data/adb must be visible inside for Debian tools)
if ! mount | grep -q '/data/local/debian/data/adb/mprime '; then
    mkdir -p /data/local/debian/data/adb/mprime
    mount -o bind /data/adb/mprime /data/local/debian/data/adb/mprime 2>/dev/null || \
        /bin/busybox mount -o bind /data/adb/mprime /data/local/debian/data/adb/mprime 2>/dev/null
fi

# 3) ddns-go (otherwise the DDNS record goes stale on next ISP re-dial)
if ! pgrep -x ddns-go >/dev/null 2>&1; then
    (cd /data/local/ddns-go && nohup /data/local/ddns-go/ddns-go >/dev/null 2>&1 &)
fi

# 4) mprime farm daemon
if ! pgrep -x mprime >/dev/null 2>&1; then
    (cd $FARM && nohup ./mprime -d >> run.log 2>&1 &)
fi

# 5) watchdog (keeps mprime/ddns-go/apn-loop alive, rotates logs)
if ! pgrep -f 'farm_watchdog.sh' >/dev/null 2>&1; then
    nohup $MP/farm_watchdog.sh >/dev/null 2>&1 &
fi

# 6) PrimeNet result submission through AutoPrimeNet (v5 API, report-only).
#    Runs on the box itself: the arm64 Debian chroot (/data/local/debian64) has
#    Python 3.11 + requests, and /data/adb/mprime is bind-mounted into it.
D64=/data/local/debian64
if [ -d "$D64" ]; then
    if ! mount | grep -q "$D64/data/adb/mprime "; then
        mkdir -p "$D64/data/adb/mprime"
        mount -o bind /data/adb/mprime "$D64/data/adb/mprime" 2>/dev/null || \
            /bin/busybox mount -o bind /data/adb/mprime "$D64/data/adb/mprime" 2>/dev/null
    fi
    if [ -x "$MP/apn/apn_loop.sh" ] && ! pgrep -f "$MP/apn/apn_loop.sh" >/dev/null 2>&1; then
        nohup $MP/apn/apn_loop.sh >/dev/null 2>&1 &
    fi
fi
exit 0
