#!/system/bin/sh
# 2 single-threaded ECM workers on the phone (the farm's real config shape)
cd /data/local/tmp/pt || exit 1
printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=2\nNumWorkers=2\nCoresPerWorker=1\nMemory=384\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' > prime.txt
printf 'ECM=1,2,264169,-1,50000,500000,1\nECM=1,2,264169,-1,50000,500000,1\n' > worktodo.txt
rm -f run.log results.txt results.json.txt ./*.bu
echo "---- 2 workers x 1 thread, ECM M264169 B1=50000, 180s"
./mprime -d > run.log 2>&1 &
mp=$!
n=0
while [ "$n" -lt 180 ]; do
  sleep 30
  n=$((n+30))
  echo "   [${n}s] lines=$(wc -l < run.log) last=$(tail -1 run.log)"
done
kill -9 "$mp" 2>/dev/null
wait "$mp" 2>/dev/null
echo "---- counters:"
for pat in "Point no longer on Edwards" "Hardware error" "Gerbicz error check" "roundoff" "FFT error" "disagrees" "stage 1 at exponent bit" "stage 2 at"; do
  echo "   [$pat] x$(grep -c "$pat" run.log 2>/dev/null)"
done
echo "---- thread lines:"
grep -E "using|threads" run.log | head -4
echo "---- last 3:"
tail -3 run.log
echo "DONE"
