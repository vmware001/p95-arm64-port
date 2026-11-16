#!/system/bin/sh
# ecmt2.sh <threads> <stage2extra> <exponent> <maxsecs>
# Tests "MT FFT for stage 1, single-threaded stage 2" via Stage2ExtraThreads.
# stage2_threads = stage1_threads + Stage2ExtraThreads, so threads=4 and
# stage2extra=-3 should give 4-thread stage 1 and 1-thread polymult (no helpers).
N=${1:-4}; EXTRA=${2:--3}; EXP=${3:-264169}; MAXS=${4:-240}
SRC=/data/local/tmp/pt/mprime
D=/data/local/tmp/ecmt2_${N}_${EXTRA}
rm -rf "$D"; mkdir -p "$D"; cp "$SRC" "$D/mprime"; chmod 755 "$D/mprime"
cd "$D" || exit 1
printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=%s\nNumWorkers=1\nCoresPerWorker=%s\nStage2ExtraThreads=%s\nMemory=512\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' "$N" "$N" "$EXTRA" > prime.txt
printf 'ECM=1,2,%s,-1,50000,500000,1\n' "$EXP" > worktodo.txt
echo "==== threads=$N Stage2ExtraThreads=$EXTRA exp=$EXP maxs=$MAXS"
./mprime -d > run.log 2>&1 &
mp=$!
n=0; prev=0; stalled=0
while [ "$n" -lt "$MAXS" ]; do
  sleep 10
  n=$((n+10))
  L=$(wc -l < run.log 2>/dev/null)
  echo "   [${n}s] lines=${L} last=$(tail -1 run.log)"
  if [ "$L" = "$prev" ]; then stalled=$((stalled+10)); else stalled=0; fi
  prev=$L
  if [ "$stalled" -ge 60 ]; then echo "   !! no log growth for 60s (hung)"; break; fi
  if grep -q "No work to do" run.log 2>/dev/null; then echo "   completed"; break; fi
done
kill -9 "$mp" 2>/dev/null
wait "$mp" 2>/dev/null
echo "---- stage2 lines : $(grep -c 'stage 2 at' run.log)"
echo "---- polymult helper lines: $(grep -c 'polymult helper' run.log)"
echo "---- ecm complete  : $(grep -c 'completed 1 ECM curve' run.log)"
echo "---- errors: edwards=$(grep -c 'Point no longer on Edwards' run.log) hw=$(grep -c 'Hardware error' run.log)"
echo "---- tail:"; tail -4 run.log
cd /; rm -rf "$D"
echo DONE
