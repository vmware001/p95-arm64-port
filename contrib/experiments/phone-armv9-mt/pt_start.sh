#!/system/bin/sh
# pt_start.sh <tag> <threads> <exponent>  -- start a detached mprime run
# Writes into /data/local/tmp/run_<tag>/ and returns immediately, so a dropped
# adb connection does not kill the test.
TAG=${1:-t1}; N=${2:-4}; EXP=${3:-264169}
D=/data/local/tmp/run_$TAG
rm -rf "$D"; mkdir -p "$D"
cp /data/local/tmp/pt/mprime "$D/mprime"; chmod 755 "$D/mprime"
cd "$D" || exit 1
printf 'UsePrimenet=0\nStressTester=0\nNumCPUs=%s\nNumWorkers=1\nCoresPerWorker=%s\nMemory=512\nProofPower=0\nV30OptionsConverted=1\nOutputIterations=0\n' "$N" "$N" > prime.txt
printf 'ECM=1,2,%s,-1,50000,500000,1\n' "$EXP" > worktodo.txt
setsid ./mprime -d > run.log 2>&1 &
sleep 2
echo "started tag=$TAG threads=$N exp=$EXP dir=$D pid=$(pgrep -f "$D/mprime" | head -1)"
