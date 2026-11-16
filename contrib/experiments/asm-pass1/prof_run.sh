#!/system/bin/sh
# PC-sample profiling run. Usage: sh prof_run.sh <bin> <exp> <iters> <usec>
B=/data/local/tmp/bench
BIN=$1
EXP=$2
ITERS=$3
USEC=$4
PID=$(cat /data/adb/mprime/farm1/mprime.pid 2>/dev/null)
echo "pid=$PID bin=$BIN exp=$EXP iters=$ITERS usec=$USEC"
[ -n "$PID" ] && kill -STOP $PID 2>/dev/null
trap 'kill -CONT $PID 2>/dev/null' EXIT INT TERM HUP
sleep 1
echo "--- governor: $(cat /sys/devices/system/cpu/cpu4/cpufreq/scaling_governor 2>/dev/null) freq: $(cat /sys/devices/system/cpu/cpu4/cpufreq/scaling_cur_freq 2>/dev/null)"
PROFILE_US=$USEC taskset 4 $B/$BIN $EXP $ITERS 2>/data/local/tmp/bench/prof_out.txt
echo "--- stdout above; sample lines: $(grep -c '^PC' /data/local/tmp/bench/prof_out.txt)"
grep -E '^(PROF|timeA|timeB|digest|profile_digest|fftlen)' /data/local/tmp/bench/prof_out.txt
kill -CONT $PID 2>/dev/null
sleep 1
echo "resumed"
