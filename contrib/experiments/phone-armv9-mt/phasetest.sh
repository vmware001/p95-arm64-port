#!/system/bin/sh
# validate heartbeat + SIGKILL phase logic on Android (box) with a local mprime
cd /data/local/tmp/pt || exit 1
printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=1\nNumWorkers=1\nCoresPerWorker=1\nMemory=128\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' > prime.txt
printf 'ECM=1,2,264169,-1,50000,500000,1\n' > worktodo.txt
rm -f run.log results.txt
./mprime -d > run.log 2>&1 &
mp=$!
n=0
while [ $n -lt 20 ]; do
  sleep 10
  n=$((n+10))
  echo "[${n}s] $(tail -1 run.log)"
done
kill -9 $mp 2>/dev/null
wait $mp 2>/dev/null
echo killed_ok
echo "progress_lines=$(grep -c 'stage 1 at exponent bit' run.log)"
echo "still_running=$(pgrep -x mprime | wc -l)"
