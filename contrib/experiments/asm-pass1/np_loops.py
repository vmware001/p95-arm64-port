import sys
lines = open('dis.txt').read().splitlines()
rows = []
for ln in lines:
    mark = ln[0]
    c = int(ln[1:6])
    rest = ln[7:]
    rows.append((mark, c, rest))
# find branches with targets in function; print loop candidates (backward branches)
addr = {}
for i, (m, c, rest) in enumerate(rows):
    a = int(rest.split()[0], 16)
    addr[a] = i
loops = []
for i, (m, c, rest) in enumerate(rows):
    p = rest.split(None, 2)
    if len(p) < 3:
        continue
    mn, ops = p[1], p[2]
    if mn.startswith('b') and '#' in ops:
        tgt = ops.split('#')[1].split()[0]
        try:
            t = int(tgt, 16)
        except ValueError:
            continue
        if t in addr and addr[t] < i:
            loops.append((addr[t], i, mn, ops))
for lo, hi, mn, ops in loops:
    body = hi - lo + 1
    samp = sum(c for m, c, r in rows[lo:hi + 1])
    a0 = int(rows[lo][2].split()[0], 16)
    a1 = int(rows[hi][2].split()[0], 16)
    print('loop %08x..%08x  instrs=%4d  samples=%3d  (%s %s)' % (a0, a1, body, samp, mn, ops))
