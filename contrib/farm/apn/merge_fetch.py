#!/usr/bin/env python3
"""Merge assignments fetched by AutoPrimeNet into the farm worktodo (worktodo.add).

AutoPrimeNet runs in a separate ``pool`` directory and writes mprime-compatible
``ECM2=`` lines there.  The farm uses one worker section per core and one curve
per worktodo line, so this script:

* splits every fetched assignment into single-curve lines (50 curves -> 50 lines,
  so the four workers can share an assignment / load-balance),
* appends them to ``worktodo.add`` grouped by worker section, always filling the
  currently least-loaded worker first,
* archives the fetched assignment lines and clears the pool.

mprime merges ``worktodo.add`` into ``worktodo.txt`` on its next save-files timer
(<= 10 minutes) and picks the new lines up without a restart.

Runs inside the arm64 Debian chroot, where /data/adb/mprime is bind-mounted.
"""
import datetime
import os
import re

MP = "/data/adb/mprime"
POOL = os.path.join(MP, "apn", "pool")
FARM = os.path.join(MP, "farm1")
WORKTODO = os.path.join(FARM, "worktodo.txt")
WORKTODO_ADD = os.path.join(FARM, "worktodo.add")
POOL_WORKTODO = os.path.join(POOL, "worktodo.txt")
ARCHIVE = os.path.join(POOL, "fetched.json.txt")

TASK_RE = re.compile(r"^(ECM2?|PRP|Pminus1|Pfactor|Test|DoubleCheck|Factor)=")
SECTION_RE = re.compile(r"^\[(Worker #\d+)\]")


def read_lines(path):
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8", errors="replace") as f:
        return [line.rstrip("\r\n") for line in f]


def section_depths(paths):
    """Return the section order and the number of task lines in each section."""
    order, counts, current = [], {}, None
    for path in paths:
        for line in read_lines(path):
            match = SECTION_RE.match(line)
            if match:
                current = match.group(1)
                if current not in counts:
                    counts[current] = 0
                    order.append(current)
            elif current and TASK_RE.match(line):
                counts[current] += 1
    return order, counts


def split_curves(line):
    """Split ECM2=aid,k,b,n,c,B1,B2,curves into single-curve lines."""
    if not line.startswith("ECM2="):
        return [line]
    parts = line.split(",")
    if len(parts) < 8:
        return [line]
    try:
        curves = int(parts[7])
    except ValueError:
        return [line]
    if curves <= 1:
        return [line]
    return [",".join(parts[:7] + ["1"]) for _ in range(curves)]


def main():
    fetched = [line.strip() for line in read_lines(POOL_WORKTODO) if TASK_RE.match(line.strip())]
    if not fetched:
        print("merge: no fetched assignments in %s" % POOL_WORKTODO)
        return 0

    new_lines = []
    for line in fetched:
        new_lines.extend(split_curves(line))
    order, counts = section_depths([WORKTODO, WORKTODO_ADD])
    if not order:
        order, counts = ["Worker #1"], {"Worker #1": 0}

    buckets = {section: [] for section in order}
    load = dict(counts)
    for line in new_lines:
        section = min(order, key=lambda name: load[name])
        buckets[section].append(line)
        load[section] += 1

    with open(WORKTODO_ADD, "a", encoding="utf-8", newline="\n") as f:
        for section in order:
            if not buckets[section]:
                continue
            f.write("[%s]\n" % section)
            for line in buckets[section]:
                f.write(line + "\n")

    stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(ARCHIVE, "a", encoding="utf-8", newline="\n") as f:
        for line in fetched:
            f.write("%s\t%s\n" % (stamp, line))
    with open(POOL_WORKTODO, "w", encoding="utf-8") as f:
        pass
    results = os.path.join(POOL, "results.json.txt")
    if os.path.exists(results):
        open(results, "w", encoding="utf-8").close()

    print(
        "merge: %d assignment(s) -> %d curve line(s) -> %s (%s)"
        % (
            len(fetched),
            len(new_lines),
            ", ".join("%s +%d" % (section, len(buckets[section])) for section in order),
            os.path.basename(WORKTODO_ADD),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
