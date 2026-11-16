#!/system/bin/sh
# Multi-run 8-thread ECM validation loop (fixed binary).
# usage: sh runfix2.sh <num_runs> <exponent> [binary]
# Result per run: ERROR (ed_check corruption), STAGE2 (stage1 completed ok), TIMEOUT.
N=${1:-2}
EXP=${2:-264169}
BIN=${3:-/data/local/tmp/mtbox/mprime-mc4mt}
B=/data/local/tmp/mtbox/ecmfix2
mkdir -p $B
cd $B || exit 1
cp -f "$BIN" ./mprime || exit 1
chmod 755 ./mprime
rm -f summary.txt
echo "exp=$EXP bin=$(md5sum ./mprime | cut -d' ' -f1)" >> summary.txt
i=1
while [ $i -le $N ]; do
  rm -f e* *.bu run.log results.txt results.json.txt mprime.pid
  printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=8\nNumWorkers=1\nCoresPerWorker=8\nMemory=256\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' > prime.txt
  printf 'ECM=1,2,%s,-1,50000,500000,1\n' "$EXP" > worktodo.txt
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
  fl=$(grep -m1 "bits-per-word below FFT limit" run.log 2>/dev/null)
  echo "RUN $i exp=$EXP result=$res elapsed=${n}s hw=$hw threads=$tid $fl" >> summary.txt
  i=$((i+1))
done
echo DONE >> summary.txt
