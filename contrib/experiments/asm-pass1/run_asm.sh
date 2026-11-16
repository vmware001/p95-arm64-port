#!/system/bin/sh
# asm-port benchmark runner.  Usage: sh run_asm.sh <tag> <bin> <exp> <iters> [rep]
# Freezes the mprime farm while measuring and ALWAYS resumes it (trap).
B=/data/local/tmp/bench
TAG=$1
BIN=$2
EXP=$3
ITERS=$4
REP=${5:-1}
PID=$(cat /data/adb/mprime/farm1/mprime.pid 2>/dev/null)
echo "mprime_pid=$PID tag=$TAG bin=$BIN exp=$EXP iters=$ITERS rep=$REP"
[ -n "$PID" ] && kill -STOP $PID 2>/dev/null
trap 'kill -CONT $PID 2>/dev/null' EXIT INT TERM HUP
sleep 1
echo -n "state_frozen: "; grep -m1 State /proc/$PID/status 2>/dev/null
echo "run_begin"
taskset 4 $B/$BIN $EXP $ITERS 2>/dev/null
echo "run_end"
kill -CONT $PID 2>/dev/null
sleep 1
echo -n "state_after: "; grep -m1 State /proc/$PID/status 2>/dev/null
echo -n "load_after: "; cat /proc/loadavg
