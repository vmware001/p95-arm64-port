#!/system/bin/sh
# mprime_known_results.sh -- end-to-end validation with SMALL EXPONENTS whose
# final primality answer is already known from the historical record.
#
# Run this on the box with the mprime binary that is to be validated. It
# creates an isolated scratch directory (no PrimeNet, no farm files touched),
# runs PRP tests on known exponents, and checks the answers:
#
#   M1277   -> COMPOSITE  (status "C")
#   M11213  -> PRIME      (status "P", a Mersenne prime found in 1963)
#   M19937  -> PRIME      (status "P", a Mersenne prime found in 1971)
#
# A correct FFT/carry/modular-reduction chain reproduces these answers. Any
# optimisation that breaks exactness will (a) produce the wrong status, or
# (b) trigger mprime's own roundoff/FFT-error detection.
#
# Usage: mprime_known_results.sh [path-to-mprime]   (default: farm binary)
set -u
BIN=${1:-/data/adb/mprime/farm1/mprime}
V=/data/local/tmp/val
rm -rf $V
mkdir -p $V || exit 1
cp "$BIN" $V/mprime || { echo "cannot copy $BIN"; exit 1; }
chmod 755 $V/mprime
cp /data/adb/mprime/farm1/local.txt $V/local.txt 2>/dev/null

cat > $V/prime.txt <<'EOF'
UsePrimenet=0
StressTester=0
NumCPUs=1
NumWorkers=1
CoresPerWorker=1
Memory=256
V30OptionsConverted=1
ProofPower=6
EOF

cat > $V/worktodo.txt <<'EOF'
PRP=1,2,1277,-1,72,0
PRP=1,2,11213,-1,72,0
PRP=1,2,19937,-1,72,0
EOF

cd $V || exit 1
timeout 300 ./mprime -d > run.log 2>&1
rc=$?

echo "=== run rc=$rc"
echo "=== results.txt"
cat results.txt 2>/dev/null
echo "=== status counts (expect 1 x C, 2 x P)"
grep -o '"status":"[A-Z]*"' results.json.txt 2>/dev/null | sort | uniq -c
echo "=== errors reported by mprime (expect none)"
grep -iE "error|roundoff|disagrees|sum(INPUTS)|fatal" run.log 2>/dev/null | head -10
echo "=== run.log tail"
tail -12 run.log 2>/dev/null
echo "=== done"
