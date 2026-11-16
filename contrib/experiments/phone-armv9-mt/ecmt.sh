#!/system/bin/sh
# ecmt.sh <threads> <exponent> <maxsecs>
# Fast reproducer for the ECM stage-2 multithread hang: a tiny exponent makes
# stage 1 finish in seconds, so the stage-2 polymult helper-thread startup is
# reached immediately (that is where mprime hung on the phone).
# Use a SMALL PRIME exponent (e.g. 1103) so mprime accepts the assignment.
N=${1:-8}; EXP=${2:-1103}; MAXS=${3:-120}
SRC=/data/local/tmp/pt/mprime
D=/data/local/tmp/ecmt_${N}_${EXP}
rm -rf "$D"; mkdir -p "$D"; cp "$SRC" "$D/mprime"; chmod 755 "$D/mprime"
cd "$D" || exit 1
printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=%s\nNumWorkers=1\nCoresPerWorker=%s\nMemory=512\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' "$N" "$N" > prime.txt
printf 'ECM=1,2,%s,-1,50000,500000,1\n' "$EXP" > worktodo.txt
echo "==== threads=$N exp=$EXP maxs=$MAXS"
./mprime -d > run.log 2>&1 &
mp=$!
n=0; stalled=0
while [ "$n" -lt "$MAXS" ]; do
  sleep 10
  n=$((n+10))
  L1=$(wc -l < run.log 2>/dev/null)
  echo "   [${n}s] lines=${L1} last=$(tail -1 run.log)"
  if [ "$L1" = "${L1_prev:-0}" ]; then stalled=$((stalled+10)); else stalled=0; fi
  L1_prev=$L1
  if [ "$stalled" -ge 60 ]; then echo "   !! no log growth for 60s (hung?)"; break; fi
  if grep -q "No work to do" run.log 2>/dev/null; then echo "   completed"; break; fi
done
kill -9 "$mp" 2>/dev/null
wait "$mp" 2>/dev/null
echo "---- stage1 lines : $(grep -c 'stage 1 at exponent bit' run.log)"
echo "---- stage2 lines : $(grep -c 'stage 2 at' run.log)"
echo "---- ecm complete : $(grep -c 'completed 1 ECM curve' run.log)"
echo "---- errors: edwards=$(grep -c 'Point no longer on Edwards' run.log) hw=$(grep -c 'Hardware error' run.log) roundoff=$(grep -c 'roundoff' run.log)"
echo "---- tail:"; tail -4 run.log
cd /; rm -rf "$D"
echo DONE
