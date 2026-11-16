#!/system/bin/sh
# Reproduce attempt on the RK3528 box: ECM with an 8-thread worker
# (the box is the same ARM64 build that failed on the Dimensity phone,
#  but until now the box had only been tested with 8-thread PRP).
cd /data/local/tmp/pt || exit 1
printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=8\nNumWorkers=1\nCoresPerWorker=8\nMemory=256\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' > prime.txt
printf 'ECM=1,2,264169,-1,50000,500000,1\n' > worktodo.txt
rm -f run.log results.txt results.json.txt ./*.bu
echo "---- ECM with 8-thread worker on 4-core box"
./mprime -d > run.log 2>&1 &
mp=$!
n=0
while [ "$n" -lt 240 ]; do
  sleep 20
  n=$((n+20))
  echo "   [${n}s] lines=$(wc -l < run.log) last=$(tail -1 run.log)"
done
kill -9 "$mp" 2>/dev/null
wait "$mp" 2>/dev/null
echo "---- counters:"
for pat in "Edwards" "Hardware error" "Gerbicz error check passed" "Gerbicz error check failed" "roundoff" "FFT error" "disagrees"; do
  echo "   [$pat] x$(cat run.log results.txt 2>/dev/null | grep -c "$pat")"
done
echo "---- threads line: $(grep -m1 -E 'threads|FFT length' run.log)"
echo "---- progress lines: $(grep -c 'stage 1 at exponent bit' run.log)"
echo "DONE"
