import sys, struct, collections

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
        secs.append(dict(name_off=name, typ=typ, addr=addr, offset=offset, size=size, link=link, entsize=entsize))
    shstr = secs[e_shstrndx]
    def sname(no):
        s = d[shstr['offset'] + no:]
        return s[:s.index(b'\0')].decode()
    for s in secs:
        s['name'] = sname(s['name_off'])
    return d, secs

def symbols(d, secs):
    syms = []
    for s in secs:
        if s['typ'] == 2:  # SYMTAB
            strtab = secs[s['link']]
            n = s['size'] // 24
            for i in range(n):
                off = s['offset'] + i * 24
                st_name, st_info, st_other, st_shndx, st_value, st_size = struct.unpack_from('<IBBHQQ', d, off)
                if st_value == 0 or st_shndx == 0:
                    continue
                nm = d[strtab['offset'] + st_name:]
                nm = nm[:nm.index(b'\0')].decode()
                syms.append((st_value, st_size, nm, st_info & 0xf))
    syms.sort()
    return syms

def main():
    binpath, profpath = sys.argv[1], sys.argv[2]
    d, secs = load_elf(binpath)
    syms = symbols(d, secs)
    starts = [s[0] for s in syms]
    import bisect
    pcs = []
    raw = open(profpath, 'rb').read()
    if raw[:2] in (b'\xff\xfe', b'\xfe\xff'):
        text = raw.decode('utf-16')
    else:
        text = raw.decode('utf-8', 'replace')
    for line in text.splitlines():
        if line.startswith('PC '):
            pcs.append(int(line.split()[1], 16))
    print('total samples', len(pcs))
    agg = collections.Counter()
    for pc in pcs:
        i = bisect.bisect_right(starts, pc) - 1
        if i < 0:
            agg['<unknown>'] += 1
            continue
        a, sz, nm, typ = syms[i]
        agg[nm] += 1
    tot = len(pcs)
    for nm, c in agg.most_common(12):
        print('%7.2f%%  %6d  %s' % (100.0 * c / tot, c, nm))
    # detailed PC histogram
    hot = collections.Counter(pcs)
    print('\ntop 40 PCs:')
    for pc, c in hot.most_common(40):
        i = bisect.bisect_right(starts, pc) - 1
        nm = syms[i][2] if i >= 0 else '?'
        print('  %08x  %5d  %s+0x%x' % (pc, c, nm, pc - syms[i][0] if i >= 0 else 0))
    # dump hot pc list per function
    with open('hot_pcs.txt', 'w') as f:
        for pc, c in sorted(hot.items()):
            i = bisect.bisect_right(starts, pc) - 1
            nm = syms[i][2] if i >= 0 else '?'
            base = syms[i][0] if i >= 0 else 0
            f.write('%x %d %s %x\n' % (pc, c, nm, pc - base if i >= 0 else 0))
    print('\nwrote hot_pcs.txt')

main()
