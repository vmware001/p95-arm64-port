#!/usr/bin/env python3
"""Rebalance the farm worktodo: keep only the selected assignments, spread their
curves evenly over the worker sections, and back the old file up.

Used when the assignment queue has to be rebuilt, e.g. because PrimeNet closed
some ECM assignments after their first curve was reported (see
contrib/primenet/patches/autoprimenet-2.0.2-ecm-keep-open.patch): their remaining
lines can no longer be credited, so they are dropped here.

Run it with mprime stopped, then start mprime again (it re-reads worktodo.txt):

    python3 rebalance_worktodo.py --keep-aid AID1 --keep-aid AID2 \
        --workers 4 --worktodo /data/adb/mprime/farm1/worktodo.txt

The script can run inside the arm64 Debian chroot (paths use /data/adb/mprime)
or directly on a machine where the file is visible.
"""
import argparse
import os
import re
import shutil
import time

TASK_RE = re.compile(r"^ECM2?=")
SECTION_RE = re.compile(r"^\[(Worker #\d+)\]")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worktodo", default="/data/adb/mprime/farm1/worktodo.txt")
    parser.add_argument("--workers", type=int, default=4, help="number of worker sections")
    parser.add_argument(
        "--keep-aid",
        action="append",
        default=[],
        metavar="AID",
        help="assignment ID to keep (repeatable); all other task lines are dropped",
    )
    args = parser.parse_args()

    if not args.keep_aid:
        raise SystemExit("give at least one --keep-aid")

    keep = set(args.keep_aid)
    workers = ["Worker #%d" % i for i in range(1, args.workers + 1)]

    with open(args.worktodo, encoding="utf-8", errors="replace") as f:
        lines = [line.rstrip("\r\n") for line in f]
    tasks = [line for line in lines if TASK_RE.match(line)]
    kept = [line for line in tasks if line.split("=", 1)[1].split(",", 1)[0] in keep]
    if not kept:
        raise SystemExit("no task lines match the given --keep-aid values")

    backup = args.worktodo + ".pre-rebalance-" + time.strftime("%Y%m%d-%H%M%S") + ".txt"
    shutil.copyfile(args.worktodo, backup)

    buckets = {w: [] for w in workers}
    for i, line in enumerate(kept):
        buckets[workers[i % len(workers)]].append(line)

    with open(args.worktodo + ".new", "w", encoding="utf-8", newline="\n") as f:
        for w in workers:
            f.write("[%s]\n" % w)
            for line in buckets[w]:
                f.write(line + "\n")
    os.replace(args.worktodo + ".new", args.worktodo)

    print("kept %d of %d task lines (%d assignment(s))" % (len(kept), len(tasks), len(keep)))
    for w in workers:
        print("  %s: %d lines" % (w, len(buckets[w])))
    print("backup:", backup)


if __name__ == "__main__":
    main()
