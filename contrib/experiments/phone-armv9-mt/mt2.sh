#!/system/bin/sh
# mt2.sh <threads> <mask> <maxsecs>
# Runs a PRP test to COMPLETION (or timeout) and reports the final error summary.
# The corruption only shows up at the final Gerbicz check / rollback, so a
# fixed 60 s slice is not enough.
N=${1:-8}; MASK=${2:-ff}; MAXS=${3:-240}
SRC=/data/local/tmp/pt/mprime
D=/data/local/tmp/mt2_${N}_${MASK}
rm -rf "$D"; mkdir -p "$D"; cp "$SRC" "$D/mprime"; chmod 755 "$D/mprime"
cd "$D" || exit 1
printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=%s\nNumWorkers=1\nCoresPerWorker=%s\nMemory=384\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=1000000\n' "$N" "$N" > prime.txt
printf 'PRP=1,2,264169,-1,72,0\n' > worktodo.txt
taskset "$MASK" ./mprime -d > run.log 2>&1 &
mp=$!
n=0; done_flag=0
while [ "$n" -lt "$MAXS" ]; do
  sleep 10
  n=$((n+10))
  if grep -q "No work to do" run.log 2>/dev/null; then done_flag=1; break; fi
done
kill -9 "$mp" 2>/dev/null
wait "$mp" 2>/dev/null
echo "==== N=$N mask=$MASK elapsed=${n}s completed=$done_flag"
echo "cfg        : $(grep -m1 'threads' run.log)"
echo "aff_fail   : $(grep -c 'Error setting affinity' run.log)"
echo "hw_err     : $(grep -c 'Hardware error' run.log)"
echo "gerb_fail  : $(grep -c 'Gerbicz error check failed' run.log)"
echo "gerb_ok    : $(grep -c 'Gerbicz error check passed' run.log)"
echo "iterations : $(grep -c 'Iteration:' run.log)"
echo "confidence : $(grep -m1 'Confidence in final result' run.log)"
echo "results    : $(tail -1 results.txt 2>/dev/null)"
cd /; rm -rf "$D"
echo DONE
