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

# 2) bind the farm dir into the debian chroot (Python lives there; the
#    linux-env wrappers chroot, so /data/adb must be visible inside)
if ! mount | grep -q '/data/local/debian/data/adb/mprime '; then
    mkdir -p /data/local/debian/data/adb/mprime
    mount -o bind /data/adb/mprime /data/local/debian/data/adb/mprime 2>/dev/null || \
        /bin/busybox mount -o bind /data/adb/mprime /data/local/debian/data/adb/mprime 2>/dev/null
fi

# 3) ddns-go (otherwise the DDNS record goes stale on next ISP re-dial)
if ! pgrep -f '/data/local/ddns-go/ddns-go' >/dev/null 2>&1; then
    (cd /data/local/ddns-go && nohup /data/local/ddns-go/ddns-go >/dev/null 2>&1 &)
fi

# 4) mprime farm daemon
if ! pgrep -x mprime >/dev/null 2>&1; then
    (cd $FARM && nohup ./mprime -d >> run.log 2>&1 &)
fi

# 5) watchdog (submits results, tops up work, keeps everything alive)
if ! pgrep -f 'farm_watchdog.sh' >/dev/null 2>&1; then
    nohup $MP/farm_watchdog.sh >/dev/null 2>&1 &
fi
exit 0
