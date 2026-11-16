import sys, struct, collections, bisect

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
        if s['typ'] == 2:
            strtab = secs[s['link']]
            n = s['size'] // 24
            for i in range(n):
                off = s['offset'] + i * 24
                st_name, st_info, st_other, st_shndx, st_value, st_size = struct.unpack_from('<IBBHQQ', d, off)
                if st_value == 0 or st_shndx == 0:
                    continue
                nm = d[strtab['offset'] + st_name:]
                nm = nm[:nm.index(b'\0')].decode()
                if st_size:
                    syms.append((st_value, st_size, nm))
    syms.sort()
    return syms

def main():
    import capstone
    binpath, func = sys.argv[1], sys.argv[2]
    pclist = sys.argv[3] if len(sys.argv) > 3 else None
    d, secs = load_elf(binpath)
    syms = symbols(d, secs)
    hit = [s for s in syms if s[2] == func]
    if not hit:
        print('symbol not found'); return
    base, size, _ = hit[0]
    print('%s at %x size %d (0x%x)' % (func, base, size, size))
    # map vaddr -> file offset using section containing base
    text = None
    for s in secs:
        if s['addr'] <= base < s['addr'] + s['size']:
            text = s
            break
    fo = text['offset'] + (base - text['addr'])
    code = d[fo:fo + size]
    md = capstone.Cs(capstone.CS_ARCH_ARM64, capstone.CS_MODE_ARM)
    md.detail = False
    insns = [(i.address, i.mnemonic, i.op_str, i.size) for i in md.disasm(code, base)]
    bname = {a: (m + ' ' + o) for a, m, o, sz in insns}
    order = [a for a, m, o, sz in insns]
    if len(sys.argv) > 4 and sys.argv[4] == 'count':
        # count instructions by mnemonic in whole function
        c = collections.Counter(m.split('.')[0] for a, m, o, sz in insns)
        tot = sum(c.values())
        print('instructions: %d' % tot)
        for k, v in c.most_common(20):
            print('  %-10s %6d  %5.1f%%' % (k, v, 100.0 * v / tot))
        return
    pcs = []
    if pclist:
        raw = open(pclist, 'rb').read()
        text2 = raw.decode('utf-16') if raw[:2] in (b'\xff\xfe', b'\xfe\xff') else raw.decode('utf-8', 'replace')
        for line in text2.splitlines():
            if line.startswith('PC '):
                pcs.append(int(line.split()[1], 16))
    hot = collections.Counter(p for p in pcs if base <= p < base + size)
    print('samples in function: %d' % sum(hot.values()))
    for pc, c in hot.most_common(25):
        idx = bisect.bisect_right(order, pc) - 1
        print('  +0x%-5x %4d  %s' % (pc - base, c, bname.get(order[idx], '?')))
        # show preceding 2 and following 2 instructions
        for k in range(max(0, idx - 2), min(len(order), idx + 3)):
            a = order[k]
            print('        %s+0x%-5x  %s' % (' ' if a != pc else '*', a - base, bname[a]))

main()
