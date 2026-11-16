import sys

lines = open('dis.txt').read().splitlines()
rows = []
for ln in lines:
    mark = ln[0]
    c = int(ln[1:6])
    rest = ln[7:]
    rows.append((mark, c, rest))
regions = []
cur = []
for idx, (mark, c, rest) in enumerate(rows):
    if c > 0:
        cur.append((idx, c, rest))
    else:
        if cur:
            regions.append(cur)
            cur = []
if cur:
    regions.append(cur)
tot = sum(c for m, c, r in rows)
print('total', tot)
for reg in regions:
    s = sum(c for i, c, r in reg)
    print('REGION insn %5d-%5d n=%3d samples=%4d (%5.1f%%)  first: %s' %
          (reg[0][0], reg[-1][0], len(reg), s, 100.0 * s / tot, reg[0][2].strip()[:70]))
