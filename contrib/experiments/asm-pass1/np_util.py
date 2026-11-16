import sys, struct, bisect

def load_elf(path):
    d = open(path, 'rb').read()
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
    return syms

syms = load_elf('bin/bench_nostrip')
starts = [s[0] for s in syms]
if sys.argv[1] == 'resolve':
    for a in [int(x, 16) for x in sys.argv[2:]]:
        i = bisect.bisect_right(starts, a) - 1
        print('%x -> %s at %x size %d off +0x%x' % (a, syms[i][2], syms[i][0], syms[i][1], a - syms[i][0]))
else:
    lo, hi = int(sys.argv[2], 16), int(sys.argv[3], 16)
    for ln in open('dis.txt'):
        rest = ln[7:]
        a = int(rest.split()[0], 16)
        if lo <= a <= hi:
            print(ln.rstrip())
