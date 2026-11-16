import sys, struct, collections, bisect, re
import capstone

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
    syms = []
    for s in secs:
        if s['typ'] == 2:
            strtab = secs[s['link']]
            n = s['size'] // 24
            for i in range(n):
                off = s['offset'] + i * 24
                st_name, st_info, st_other, st_shndx, st_value, st_size = struct.unpack_from('<IBBHQQ', d, off)
                nm = d[strtab['offset'] + st_name:]
                nm = nm[:nm.index(b'\0')].decode()
                if st_value and st_size and st_shndx:
                    syms.append((st_value, st_size, nm))
    syms.sort()
    return d, secs, syms

def fn_disasm(binpath, func):
    d, secs, syms = load_elf(binpath)
    ent = [s for s in syms if s[2] == func][0]
    base, size, _ = ent
    text = None
    for s in secs:
        if s['addr'] <= base < s['addr'] + s['size']:
            text = s; break
    fo = text['offset'] + (base - text['addr'])
    code = d[fo:fo + size]
    md = capstone.Cs(capstone.CS_ARCH_ARM64, capstone.CS_MODE_ARM)
    return base, size, [(i.address, i.mnemonic, i.op_str) for i in md.disasm(code, base)]

def main():
    binpath, func, profpath = sys.argv[1], sys.argv[2], sys.argv[3]
    base, size, ins = fn_disasm(binpath, func)
    addrs = [a for a, m, o in ins]
    cnt = collections.Counter()
    raw = open(profpath, 'rb').read()
    t = raw.decode('utf-16') if raw[:2] in (b'\xff\xfe', b'\xfe\xff') else raw.decode('utf-8', 'replace')
    pcs = [int(l.split()[1], 16) for l in t.splitlines() if l.startswith('PC ')]
    for pc in pcs:
        if base <= pc < base + size:
            i = bisect.bisect_right(addrs, pc) - 1
            cnt[i] += 1
    tot = sum(cnt.values())
    # find loop back-edges: branch targets < current address
    targets = set()
    for a, m, o in ins:
        if m.startswith('b') and o.startswith('#0x'):
            try:
                targets.add(int(o[1:], 16))
            except ValueError:
                pass
    out = open('dis.txt', 'w')
    cum = 0
    for i, (a, m, o) in enumerate(ins):
        c = cnt.get(i, 0)
        cum += c
        mark = '*' if c else ' '
        tgt = '   <-- target' if a in targets else ''
        out.write('%s%5d %08x %-6s %-34s %s\n' % (mark, c, a, m, o, tgt))
    out.close()
    print('function %s base=%x size=%d instrs=%d samples=%d' % (func, base, size, len(ins), tot))

main()
