#!/system/bin/sh
# Multi-run 8-thread ECM reproduction loop for the MT corruption fix.
# usage: sh runfix.sh [num_runs]
# Each run: fresh ECM on M264169 B1=50000 with CoresPerWorker=8.
# Result: ERROR (ed_check failed = corruption), STAGE2 (stage 1 completed and
# stage 2 started = success), or TIMEOUT.
N=${1:-3}
B=/data/local/tmp/mtbox/ecmfix
mkdir -p $B
cd $B || exit 1
cp -f /data/local/tmp/mtbox/mprime-mtfix1 ./mprime || exit 1
chmod 755 ./mprime
rm -f summary.txt
echo "fix test: $(md5sum ./mprime)" >> summary.txt
i=1
while [ $i -le $N ]; do
  rm -f e* *.bu run.log results.txt results.json.txt mprime.pid
  printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=8\nNumWorkers=1\nCoresPerWorker=8\nMemory=256\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' > prime.txt
  printf 'ECM=1,2,264169,-1,50000,500000,1\n' > worktodo.txt
  (nohup ./mprime -d > run.log 2>&1 &)
  n=0
  res=TIMEOUT
  while [ $n -lt 540 ]; do
    sleep 5
    n=$((n+5))
    if grep -q "Point no longer" run.log 2>/dev/null; then res=ERROR; break; fi
    if grep -q "polymult helper thread" run.log 2>/dev/null; then res=STAGE2; break; fi
  done
  pkill -9 -x mprime 2>/dev/null
  sleep 2
  hw=$(grep -c "Hardware error" run.log 2>/dev/null)
  tid=$(grep -m1 -o "[0-9] threads" run.log 2>/dev/null)
  echo "RUN $i result=$res elapsed=${n}s hw=$hw threads=$tid" >> summary.txt
  i=$((i+1))
done
echo DONE >> summary.txt
