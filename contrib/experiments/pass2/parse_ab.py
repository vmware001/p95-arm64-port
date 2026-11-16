#!/usr/bin/env python
# parse_ab.py -- summarise ab2.sh interleaved benchmark logs.
# usage: python parse_ab.py <log> [tolerance_pct]
import sys, re, statistics, collections

log = sys.argv[1]
tol = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
rows = collections.defaultdict(lambda: collections.defaultdict(list))
dig = collections.defaultdict(set)
txt = open(log, 'r', errors='replace').read()
for m in re.finditer(r'REP (\d+) RUN (\d+) BIN (\S+) (timeA_ms|timeB_ms)=([\d.]+)', txt):
    r, i, spec, k, v = m.group(1), m.group(2), m.group(3), m.group(4), float(m.group(5))
    rows[spec][k].append(v)
for m in re.finditer(r'REP (\d+) RUN (\d+) BIN (\S+) digest=([0-9a-f]+)', txt):
    dig[m.group(3)].add(m.group(4))
fails = re.findall(r'.*FAIL.*', txt)
done = '=== done' in txt

print('log %s  finished=%s  fails=%d' % (log, done, len(fails)))
for f in fails:
    print('  ' + f.strip())

specs = list(rows.keys())
base = None
for s in specs:
    if s == 'bench_base':
        base = s
print()
hdr = '%-16s %4s %10s %10s %10s   %s' % ('spec', 'n', 'median', 'min', 'max', 'digest')
print(hdr)
print('-' * len(hdr))
for s in specs:
    for k in ('timeA_ms', 'timeB_ms'):
        v = sorted(rows[s][k])
        if not v:
            continue
        med = statistics.median(v)
        mark = ''
        if base and s != base and k == 'timeA_ms':
            bv = rows[base][k]
            if bv:
                bmed = statistics.median(bv)
                d = 100.0 * (med - bmed) / bmed
                mark = '  vs base %+0.2f%%%s' % (d, '  <-- beyond noise' if abs(d) > tol else '')
        print('%-16s %4d %10.4f %10.4f %10.4f   %s%s'
              % (s if k == 'timeA_ms' else '  ' + k, len(v), med, v[0], v[-1],
                 ','.join(sorted(dig.get(s, []))), mark))
print()
if base and base in rows:
    bmed = statistics.median(rows[base]['timeA_ms'])
    print('reference: %s median timeA = %.4f ms (%d runs)' % (base, bmed, len(rows[base]['timeA_ms'])))
