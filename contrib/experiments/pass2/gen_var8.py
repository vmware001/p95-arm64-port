#!/usr/bin/env python
# gen_var8.py -- generate restructured copies of gwnum/r4/r4pass2_8.c.
# usage: python gen_var8.py <mode> <outfile>
#   mode = loopc : unroll the 4-iteration inner butterfly loops by 2
#          big   : unroll the 32/64-iteration butterfly loops by 2
#          all   : both of the above
#          pf    : leave loops alone; enable twiddle prefetch (handled by -D flags)
# Only the complex pass-2 body is touched; per-element arithmetic is unchanged
# (each macro call is emitted in the same order with the same arguments).
import sys, re

src = 'gwnum/r4/r4pass2_8.c'
mode = sys.argv[1]
out = sys.argv[2]
lines = open(src, 'r', encoding='utf-8', errors='replace').read().split('\n')

# locate the complex pass-2 function body
start = None
for i, l in enumerate(lines):
    if 'void r4_pass2_8_levels_complex(struct gwasm_data *__restrict g) {' in l:
        start = i
        break
assert start is not None, 'complex function not found'

# (loop-pattern, macro-anchor expected within the next 4 lines, body lines,
#  new trip count, category)
EDITS = [
    ('for(unsigned int loopA = 64;',  'r4_f2cl_four_complex_djbfft',      3, 32, 'big'),
    ('for(unsigned int loopC = 4;',   'r4_x4cl_four_complex_djbfft(rsi, 4*64', 3,  2, 'loopc'),
    ('for(unsigned int loopC = 4;',   'r4_x4cl_four_complex_djbfft(rsi, 16*64', 1, 2, 'loopc'),
    ('for(unsigned int loopA = 32;',  'r4_x4cl_four_complex_with_square', 1, 16, 'big'),
    ('for(unsigned int loopC = 4;',   'r4_x4cl_four_complex_djbunfft(rsi, 8*64', 1, 2, 'loopc'),
    ('for(unsigned int loopC = 4;',   'r4_x4cl_four_complex_djbunfft(rsi, 2*64', 2, 2, 'loopc'),
    ('for(unsigned int loopA = 32;',  'r4_x4cl_four_complex_djbunfft(rsi, 64, 32*64', 1, 16, 'big'),
]

want = {'loopc': {'loopc'}, 'big': {'big'}, 'all': {'loopc', 'big'}, 'pf': set()}[mode]

pos = start
for pat, anchor, nbody, newtrip, cat in EDITS:
    # find the next loop header matching pat whose body starts with anchor
    i = pos
    found = None
    while i < len(lines):
        if pat in lines[i]:
            window = '\n'.join(lines[i + 1:i + 5])
            if anchor in window:
                found = i
                break
        i += 1
    assert found is not None, 'edit not found: %s / %s' % (pat, anchor)
    pos = found + nbody + 1
    if cat not in want:
        continue
    header = lines[found]
    oldtrip = int(re.search(r'loop[ABC] = (\d+);', header).group(1))
    assert oldtrip == newtrip * 2, (oldtrip, newtrip)
    lines[found] = re.sub(r'(loop[ABC] = )\d+(;)', r'\g<1>%d\g<2>' % newtrip, header)
    body = lines[found + 1:found + 1 + nbody]
    lines[found + 1 + nbody:found + 1 + nbody] = body   # duplicate the body

open(out, 'w', encoding='utf-8').write('\n'.join(lines))
print('wrote %s (mode=%s)' % (out, mode))
