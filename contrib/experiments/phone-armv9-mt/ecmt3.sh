#!/system/bin/sh
# ecmt3.sh <threads> <exponent> <maxsecs>
# Distinguishes "silent stage 2" from "deadlock" by sampling the process CPU
# time (/proc/<pid>/stat utime+stime) and RSS while mprime runs one full ECM
# curve. Prints a heartbeat every 20 s.
N=${1:-4}; EXP=${2:-264169}; MAXS=${3:-900}
SRC=/data/local/tmp/pt/mprime
D=/data/local/tmp/ecmt3_${N}
rm -rf "$D"; mkdir -p "$D"; cp "$SRC" "$D/mprime"; chmod 755 "$D/mprime"
cd "$D" || exit 1
printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=%s\nNumWorkers=1\nCoresPerWorker=%s\nMemory=512\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' "$N" "$N" > prime.txt
printf 'ECM=1,2,%s,-1,50000,500000,1\n' "$EXP" > worktodo.txt
echo "==== threads=$N exp=$EXP maxs=$MAXS"
./mprime -d > run.log 2>&1 &
mp=$!
cpu0=$(sed 's/.*) //' /proc/$mp/stat 2>/dev/null | awk '{print $12+$13}')
n=0; prev=0
while [ "$n" -lt "$MAXS" ]; do
  sleep 20
  n=$((n+20))
  cpu=$(sed 's/.*) //' /proc/$mp/stat 2>/dev/null | awk '{print $12+$13}')
  rss=$(awk '/VmRSS/{print $2}' /proc/$mp/status 2>/dev/null)
  dcpu=$((cpu - cpu0)); prev=$cpu
  pct=$((dcpu * 100 / 2))   # ticks(100Hz) over 20s => /(20*100)*100
  echo "   [${n}s] lines=$(wc -l < run.log) cpu=${pct}% rss=${rss}kB last=$(tail -1 run.log)"
  if grep -qE "completed 1 ECM curve|has a factor|No work to do" run.log 2>/dev/null; then echo "   FINISHED"; break; fi
done
kill -9 "$mp" 2>/dev/null
wait "$mp" 2>/dev/null
echo "---- stage1 lines : $(grep -c 'stage 1 at exponent bit' run.log)"
echo "---- stage2 lines : $(grep -c 'stage 2 at' run.log)"
echo "---- completed    : $(grep -cE 'completed 1 ECM curve|has a factor' run.log)"
echo "---- errors: edwards=$(grep -c 'Point no longer on Edwards' run.log) hw=$(grep -c 'Hardware error' run.log) roundoff=$(grep -c 'roundoff' run.log)"
echo "---- results:"; cat results.txt 2>/dev/null | tail -3
echo "---- tail:"; tail -5 run.log
cd /; rm -rf "$D"
echo DONE
