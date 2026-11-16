#!/system/bin/sh
# run.sh <phase>  — one validation phase on Android (MuMu or a real phone, root or shell)
# phases: prp1 | ecm1 | prp8 | ecm8 | all
cd /data/local/tmp/pt || exit 1

write_cfg() {
  case "$1" in
    prp1) printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=1\nNumWorkers=1\nCoresPerWorker=1\nMemory=128\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' > prime.txt
          printf 'PRP=1,2,11213,-1,72,0\nPRP=1,2,1277,-1,72,0\n' > worktodo.txt ;;
    ecm1) printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=1\nNumWorkers=1\nCoresPerWorker=1\nMemory=128\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' > prime.txt
          printf 'ECM=1,2,264169,-1,50000,500000,1\n' > worktodo.txt ;;
    prp8) printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=8\nNumWorkers=1\nCoresPerWorker=8\nMemory=256\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' > prime.txt
          printf 'PRP=1,2,264169,-1,72,0\n' > worktodo.txt ;;
    ecm8) printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=8\nNumWorkers=1\nCoresPerWorker=8\nMemory=256\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' > prime.txt
          printf 'ECM=1,2,264169,-1,50000,500000,1\n' > worktodo.txt ;;
  esac
}

phase() {
  name=$1; secs=$2; pin=$3
  write_cfg "$name"
  rm -f run.log results.txt results.json.txt ./*.bu
  echo "---- phase $name (${secs}s, pin=$pin)"
  if [ "$pin" = "1" ]; then
    # Android toybox taskset takes a bitmask: 0x80 = CPU7 (biggest core numbering)
    NI=$(nproc 2>/dev/null); [ -z "$NI" ] && NI=8
    MASK=$(printf '%x' $((1 << (NI - 1))))
    taskset "$MASK" ./mprime -d > run.log 2>&1 &
  else
    ./mprime -d > run.log 2>&1 &
  fi
  mp=$!
  n=0
  while [ "$n" -lt "$secs" ]; do
    sleep 10
    n=$((n+10))
    echo "   [${n}s] lines=$(wc -l < run.log) last=$(tail -1 run.log)"
  done
  kill -9 "$mp" 2>/dev/null
  wait "$mp" 2>/dev/null
  echo "   counters:"
  for pat in "Point no longer on Edwards" "Hardware error" "Gerbicz error check passed" "Gerbicz error check failed" "roundoff" "FFT error" "disagrees" "is a probable prime" "is not prime"; do
    echo "     [$pat] x$(cat run.log results.txt 2>/dev/null | grep -c "$pat")"
  done
  echo "   threads line: $(grep -m1 -E 'threads|FFT length' run.log)"
  echo "   results: $(tail -2 results.txt 2>/dev/null | tr '\n' ' | ')"
  cp -f run.log "run-$name.log"
}

case "${1:-all}" in
  prp1) phase prp1 40 1 ;;
  ecm1) phase ecm1 120 1 ;;
  prp8) phase prp8 120 0 ;;
  ecm8) phase ecm8 120 0 ;;
  all)  phase prp1 40 1; phase ecm1 120 1; phase prp8 120 0; phase ecm8 120 0 ;;
  *) echo "usage: sh run.sh prp1|ecm1|prp8|ecm8|all"; exit 2 ;;
esac
echo "DONE"
