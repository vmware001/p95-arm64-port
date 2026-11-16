#!/system/bin/sh
# Interleaved A/B benchmark, many short reps, farm frozen per run.
# Usage: sh ab2.sh <exp> <iters> <reps>
B=/data/local/tmp/bench
EXP=${1:-264169}
ITERS=${2:-300}
REPS=${3:-12}
PID=$(cat /data/adb/mprime/farm1/mprime.pid 2>/dev/null)
trap 'kill -CONT $PID 2>/dev/null' EXIT INT TERM HUP
echo "pid=$PID exp=$EXP iters=$ITERS reps=$REPS"
for r in $(seq 1 $REPS); do
  for v in base asm; do
    kill -STOP $PID 2>/dev/null
    sleep 0.3
    out=$(taskset 4 $B/bench_$v $EXP $ITERS 2>/dev/null)
    ta=$(echo "$out" | sed -n 's/^timeA_ms=//p')
    tb=$(echo "$out" | sed -n 's/^timeB_ms=//p')
    dg=$(echo "$out" | sed -n 's/^digest=//p')
    echo "rep$r $v A=$ta B=$tb d=$dg"
    kill -CONT $PID 2>/dev/null
  done
done
echo "=== done"
