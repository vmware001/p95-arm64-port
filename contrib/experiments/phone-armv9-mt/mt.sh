#!/system/bin/sh
# mt.sh <threads> <cpu-mask> <seconds>  -- multithread stress test on Android
# mask is a taskset bitmask (ff = all 8 cores, 0f = cpu0-3, 70 = cpu4-6, f0 = cpu4-7)
N=${1:-8}; MASK=${2:-ff}; SECS=${3:-60}
SRC=/data/local/tmp/pt/mprime
D=/data/local/tmp/mt_${N}_${MASK}
rm -rf "$D"; mkdir -p "$D"; cp "$SRC" "$D/mprime"; chmod 755 "$D/mprime"
cd "$D" || exit 1
printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=%s\nNumWorkers=1\nCoresPerWorker=%s\nMemory=384\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' "$N" "$N" > prime.txt
printf 'PRP=1,2,264169,-1,72,0\n' > worktodo.txt
echo "==== threads=$N mask=$MASK secs=$SECS"
taskset "$MASK" ./mprime -d > run.log 2>&1 &
mp=$!
sleep "$SECS"
kill -9 "$mp" 2>/dev/null
wait "$mp" 2>/dev/null
echo "cfg      : $(grep -m1 -E 'threads' run.log)"
echo "affinity_fail: $(grep -c 'Error setting affinity' run.log)"
echo "hardware_err : $(grep -c 'Hardware error' run.log)"
echo "gerbicz_fail : $(grep -c 'Gerbicz error check failed' run.log)"
echo "gerbicz_ok   : $(grep -c 'Gerbicz error check passed' run.log)"
echo "iterations   : $(grep -c 'Iteration:' run.log)"
echo "last         : $(tail -1 run.log)"
cd /; rm -rf "$D"
echo DONE
