#!/system/bin/sh
# pt_status.sh <tag> -- report progress of a detached run
TAG=${1:-t1}
D=/data/local/tmp/run_$TAG
[ -d "$D" ] || { echo "no such run: $D"; exit 1; }
pid=$(pgrep -f "$D/mprime" | head -1)
echo "dir=$D pid=${pid:-none}"
echo "lines=$(wc -l < "$D/run.log" 2>/dev/null)"
if [ -n "$pid" ]; then
  cpu=$(sed 's/.*) //' /proc/$pid/stat 2>/dev/null | awk '{print $12+$13}')
  rss=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null)
  thr=$(ls /proc/$pid/task 2>/dev/null | wc -l)
  echo "cpu_ticks=$cpu rss=${rss}kB threads=${thr}"
fi
echo "stage1=$(grep -c 'stage 1 at exponent bit' "$D/run.log" 2>/dev/null)"
echo "stage2=$(grep -c 'stage 2 at' "$D/run.log" 2>/dev/null)"
echo "completed=$(grep -cE 'completed 1 ECM curve|has a factor' "$D/run.log" 2>/dev/null)"
echo "errors: edwards=$(grep -c 'Point no longer on Edwards' "$D/run.log" 2>/dev/null) hw=$(grep -c 'Hardware error' "$D/run.log" 2>/dev/null) roundoff=$(grep -c 'roundoff' "$D/run.log" 2>/dev/null)"
echo "last: $(tail -1 "$D/run.log" 2>/dev/null)"
