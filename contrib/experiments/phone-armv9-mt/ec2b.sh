#!/system/bin/sh
# clean 2-worker single-thread ECM test in a fresh directory (no stale save files)
D=/data/local/tmp/pt2
rm -rf "$D"
mkdir -p "$D"
cp /data/local/tmp/pt/mprime "$D/mprime"
chmod 755 "$D/mprime"
cd "$D" || exit 1
printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=2\nNumWorkers=2\nCoresPerWorker=1\nMemory=512\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' > prime.txt
printf 'ECM=1,2,264169,-1,50000,500000,1\nECM=1,2,279863,-1,50000,500000,1\n' > worktodo.txt
echo "---- 2 workers x 1 thread, ECM M264169 + M279863, B1=50000, 600s"
./mprime -d > run.log 2>&1 &
mp=$!
n=0
while [ "$n" -lt 600 ]; do
  sleep 60
  n=$((n+60))
  echo "   [${n}s] lines=$(wc -l < run.log) last=$(tail -1 run.log)"
done
kill -9 "$mp" 2>/dev/null
wait "$mp" 2>/dev/null
echo "---- counters:"
for pat in "Point no longer on Edwards" "Hardware error" "roundoff" "FFT error" "disagrees" "Gerbicz error check"; do
  echo "   [$pat] x$(grep -c "$pat" run.log 2>/dev/null)"
done
echo "---- workers:"
grep -c "Worker #1" run.log
grep -c "Worker #2" run.log
echo "---- per-worker last line:"
grep "Worker #1" run.log | tail -1
grep "Worker #2" run.log | tail -1
echo "---- stage1 progress lines per worker:"
echo "   W1=$(grep 'Worker #1' run.log | grep -c 'stage 1 at exponent bit')"
echo "   W2=$(grep 'Worker #2' run.log | grep -c 'stage 1 at exponent bit')"
echo "---- stage2 lines:"
grep -c "stage 2" run.log
echo "DONE"
