#!/system/bin/sh
# Run one command pinned to cpu2 with the mprime farm frozen; always resumes.
# Usage: sh frozen.sh <cmd> [args...]
PID=$(cat /data/adb/mprime/farm1/mprime.pid 2>/dev/null)
trap 'kill -CONT $PID 2>/dev/null' EXIT INT TERM HUP
echo "pid=$PID"
kill -STOP $PID 2>/dev/null
sleep 1
timeout 120 taskset 4 "$@"
rc=$?
kill -CONT $PID 2>/dev/null
echo "rc=$rc"
