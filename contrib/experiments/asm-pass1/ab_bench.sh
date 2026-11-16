#!/system/bin/sh
# Interleaved A/B benchmark for the asm-port experiment (16K FFT).
# Usage: sh ab_bench.sh <exp> <iters> <reps>
# Every run freezes the farm and always resumes it.
B=/data/local/tmp/bench
EXP=${1:-264169}
ITERS=${2:-1000}
REPS=${3:-3}
PID=$(cat /data/adb/mprime/farm1/mprime.pid 2>/dev/null)
trap 'kill -CONT $PID 2>/dev/null' EXIT INT TERM HUP
echo "pid=$PID exp=$EXP iters=$ITERS reps=$REPS"
for r in $(seq 1 $REPS); do
  for v in base asm; do
    kill -STOP $PID 2>/dev/null
    sleep 1
    echo "--- rep$r $v"
    taskset 4 $B/bench_$v $EXP $ITERS 2>/dev/null
    kill -CONT $PID 2>/dev/null
    echo -n "load: "; cat /proc/loadavg
  done
done
echo "=== done"
