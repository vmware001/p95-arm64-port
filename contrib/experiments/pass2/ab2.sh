#!/system/bin/sh
# Interleaved A/B benchmark for the pass2 experiment.
# Usage: sh ab2.sh <exp> <iters> <nreps> <runs_per_rep> <spec...>
#   spec = binary   or   binary:variant   (binary:variant runs with GWNUM_P2VAR=variant)
# Freezes the mprime farm for the whole batch and ALWAYS resumes it.
# A hard watchdog resumes the farm after MAXSEC no matter what happens to us.
B=/data/local/tmp/bench
EXP=${1:-264169}
ITERS=${2:-1000}
REPS=${3:-3}
RUNS=${4:-10}
MAXSEC=420
shift 4
SPECS="$*"
PID=$(cat /data/adb/mprime/farm1/mprime.pid 2>/dev/null)
trap 'kill -CONT $PID 2>/dev/null' EXIT INT TERM HUP
echo "pid=$PID exp=$EXP iters=$ITERS reps=$REPS runs=$RUNS maxsec=$MAXSEC specs=$SPECS"
kill -STOP $PID 2>/dev/null
# hard watchdog: resume the farm even if this script is killed outright
( sleep $MAXSEC; kill -CONT $PID 2>/dev/null ) &
WD=$!
sleep 1
echo -n "frozen_state: "; grep -m1 State /proc/$PID/status 2>/dev/null
echo -n "freq_before: "; cat /sys/devices/system/cpu/cpu2/cpufreq/scaling_cur_freq 2>/dev/null
for r in $(seq 1 $REPS); do
  for i in $(seq 1 $RUNS); do
    for s in $SPECS; do
      b=${s%%:*}
      case "$s" in
        *:*) v=${s#*:}
             timeout 120 env GWNUM_P2VAR=$v taskset 4 $B/$b $EXP $ITERS 2>/dev/null | sed "s/^/REP $r RUN $i BIN $s /"
             rc=$? ;;
        *)   timeout 120 taskset 4 $B/$b $EXP $ITERS 2>/dev/null | sed "s/^/REP $r RUN $i BIN $s /"
             rc=$? ;;
      esac
      [ $rc -ne 0 ] && echo "REP $r RUN $i BIN $s FAIL rc=$rc"
    done
  done
done
kill $WD 2>/dev/null
echo -n "freq_after: "; cat /sys/devices/system/cpu/cpu2/cpufreq/scaling_cur_freq 2>/dev/null
echo -n "load: "; cat /proc/loadavg
kill -CONT $PID 2>/dev/null
sleep 1
echo -n "state_after: "; grep -m1 State /proc/$PID/status 2>/dev/null
echo "=== done"
