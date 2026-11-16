#!/usr/bin/env python3
"""Trim the farm results file after an AutoPrimeNet report-only run.

Lines that appear in AutoPrimeNet's results_sent.txt have been handled by the
server (credited, or definitively rejected such as duplicate "error 47"), so
move them to results.submitted.json.txt and keep everything else for the next
attempt. Runs inside the arm64 Debian chroot, where /data/adb/mprime is
bind-mounted.
"""
import os

FARM = "/data/adb/mprime/farm1"
APN = "/data/adb/mprime/apn"
RESULTS = os.path.join(FARM, "results.json.txt")
SENT = os.path.join(APN, "results_sent.txt")
ARCHIVE = os.path.join(FARM, "results.submitted.json.txt")


def read_lines(path):
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8", errors="replace") as f:
        return [line.rstrip("\r\n") for line in f if line.strip()]


def write_lines(path, lines):
    text = "\n".join(lines) + ("\n" if lines else "")
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as f:
        f.write(text)
    os.replace(tmp, path)


def main():
    sent = set(read_lines(SENT))
    current = read_lines(RESULTS)
    handled = [line for line in current if line in sent]
    remaining = [line for line in current if line not in sent]
    if handled:
        with open(ARCHIVE, "a", encoding="utf-8") as f:
            f.write("\n".join(handled) + "\n")
        write_lines(RESULTS, remaining)
    print("trim: handled=%d remaining=%d (archive=%s)" % (len(handled), len(remaining), ARCHIVE))


if __name__ == "__main__":
    main()
