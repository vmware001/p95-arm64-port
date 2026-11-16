#!/usr/bin/env python
# p2_analyze.py -- attribute SIGPROF PC samples to gwnum functions and to loops
# inside a chosen function, and report mnemonic mix per loop.
#
# usage: python p2_analyze.py <binary-with-symbols> <profile.txt> [func ...]
import sys, struct, bisect, collections, re

def load_elf(path):
    d = open(path, 'rb').read()
    assert d[:4] == b'\x7fELF'
    e_shoff = struct.unpack_from('<Q', d, 0x28)[0]
    e_shentsize = struct.unpack_from('<H', d, 0x3a)[0]
    e_shnum = struct.unpack_from('<H', d, 0x3c)[0]
    e_shstrndx = struct.unpack_from('<H', d, 0x3e)[0]
    secs = []
    for i in range(e_shnum):
        off = e_shoff + i * e_shentsize
        name, typ, flags, addr, offset, size, link, info, align, entsize = struct.unpack_from('<IIQQQQIIQQ', d, off)
        secs.append(dict(name_off=name, typ=typ, addr=addr, offset=offset, size=size, link=link))
    shstr = secs[e_shstrndx]
    def sname(no):
        s = d[shstr['offset'] + no:]
        return s[:s.index(b'\0')].decode()
    for s in secs:
        s['name'] = sname(s['name_off'])
    return d, secs

def load_syms(d, secs):
    syms = []
    for s in secs:
        if s['typ'] == 2:
            strtab = secs[s['link']]
            n = s['size'] // 24
            for i in range(n):
                off = s['offset'] + i * 24
                st_name, st_info, st_other, st_shndx, st_value, st_size = struct.unpack_from('<IBBHQQ', d, off)
                if st_value == 0 or st_shndx == 0 or st_size == 0:
                    continue
                nm = d[strtab['offset'] + st_name:]
                nm = nm[:nm.index(b'\0')].decode()
                syms.append((st_value, st_size, nm, st_info & 0xf))
    syms.sort()
    return syms

def read_profile(path):
    raw = open(path, 'rb').read()
    if raw[:2] in (b'\xff\xfe', b'\xfe\xff'):
        t = raw.decode('utf-16')
    else:
        t = raw.decode('utf-8', 'replace')
    return [int(l.split()[1], 16) for l in t.splitlines() if l.startswith('PC ')]

FP = {'fadd', 'fsub', 'fmul', 'fmla', 'fmls', 'fdiv', 'fneg', 'fabs', 'fcmp',
      'fcvt', 'fcvtzs', 'fcvtms', 'frintx', 'frintn', 'frintz', 'scvtf'}
MEM = {'ldr', 'str', 'ldp', 'stp', 'ldur', 'stur', 'ld1', 'st1', 'ldpsw'}
BR = {'b', 'bl', 'br', 'blr', 'ret', 'cbz', 'cbnz', 'tbz', 'tbnz'}

def classify(mn):
    if mn in FP: return 'fp'
    if mn in MEM: return 'mem'
    if mn.startswith('b'): return 'br'
    return 'other'

def main():
    binpath, profpath = sys.argv[1], sys.argv[2]
    funcs = sys.argv[3:]
    d, secs = load_elf(binpath)
    syms = load_syms(d, secs)
    starts = [s[0] for s in syms]
    pcs = read_profile(profpath)
    print('total samples %d' % len(pcs))

    # ---- function level attribution ----
    cnt = collections.Counter()
    nomap = 0
    for pc in pcs:
        i = bisect.bisect_right(starts, pc) - 1
        if i < 0 or pc >= syms[i][0] + syms[i][1]:
            nomap += 1
            continue
        cnt[syms[i][2]] += 1
    tot = len(pcs) - nomap
    print('unmapped %d (%.2f%%)' % (nomap, 100.0 * nomap / max(1, len(pcs))))
    for nm, c in cnt.most_common(25):
        print('  %-45s %6d  %5.1f%%' % (nm[:45], c, 100.0 * c / tot))

    # ---- loop level analysis of selected functions ----
    try:
        import capstone
    except ImportError:
        print('no capstone'); return
    md = capstone.Cs(capstone.CS_ARCH_ARM64, capstone.CS_MODE_ARM)
    for fn in funcs:
        ent = [s for s in syms if s[2] == fn]
        if not ent:
            print('function %s not found' % fn); continue
        base, size, _, _ = ent[0]
        text = None
        for s in secs:
            if s['addr'] <= base < s['addr'] + s['size']:
                text = s; break
        fo = text['offset'] + (base - text['addr'])
        code = d[fo:fo + size]
        ins = [(i.address, i.mnemonic, i.op_str) for i in md.disasm(code, base)]
        addrs = [a for a, m, o in ins]
        per = collections.Counter()
        fn_samp = 0
        for pc in pcs:
            if base <= pc < base + size:
                j = bisect.bisect_right(addrs, pc) - 1
                per[j] += 1
                fn_samp += 1
        print('\n=== %s base=%x size=%d instrs=%d samples=%d (%.1f%% of all)' %
              (fn, base, size, len(ins), fn_samp, 100.0 * fn_samp / max(1, len(pcs))))
        # back-edge loops
        idx = {a: i for i, (a, m, o) in enumerate(ins)}
        loops = []
        for i, (a, m, o) in enumerate(ins):
            if m.startswith('b') and o.startswith('#0x'):
                try: tgt = int(o[1:], 16)
                except ValueError: continue
                if tgt in idx and idx[tgt] < i:
                    loops.append((idx[tgt], i))
        # keep only innermost loops (no other loop fully inside)
        loops2 = []
        for lo, hi in loops:
            inner = [x for x in loops if x is not (lo, hi) and lo <= x[0] and x[1] <= hi and (x[0], x[1]) != (lo, hi)]
            if not inner:
                loops2.append((lo, hi))
        for lo, hi in sorted(loops2, key=lambda p: -sum(per[k] for k in range(p[0], p[1] + 1))):
            s = sum(per[k] for k in range(lo, hi + 1))
            mh = collections.Counter(classify(ins[k][1]) for k in range(lo, hi + 1))
            mns = collections.Counter(ins[k][1] for k in range(lo, hi + 1))
            print('  loop %08x..%08x body=%3d samples=%5d (%5.1f%% of fn) fp=%d mem=%d br=%d oth=%d  %s'
                  % (ins[lo][0], ins[hi][0], hi - lo + 1, s, 100.0 * s / max(1, fn_samp),
                     mh['fp'], mh['mem'], mh['br'], mh['other'],
                     ' '.join('%s:%d' % (k, v) for k, v in mns.most_common(8))))
        # residual (outside loops)
        inloop = set()
        for lo, hi in loops:
            inloop.update(range(lo, hi + 1))
        resid = sum(v for k, v in per.items() if k not in inloop)
        print('  outside-innermost-loops: samples=%d (%.1f%%)' % (resid, 100.0 * resid / max(1, fn_samp)))

main()
