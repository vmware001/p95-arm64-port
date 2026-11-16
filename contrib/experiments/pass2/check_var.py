#!/usr/bin/env python
# check_var.py -- sanity check the generated variants and the harness binary.
import sys, struct, difflib, collections
sys.path.insert(0, '.')

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
                if st_value and st_size and st_shndx:
                    nm = d[strtab['offset'] + st_name:]
                    syms.append((st_value, st_size, nm[:nm.index(b'\0')].decode()))
    syms.sort()
    return d, secs, syms

mode = sys.argv[1]
if mode == 'diff':
    name = sys.argv[2]
    a = open('gwnum/r4/r4pass2_8.c').read().split('\n')
    b = open('gwnum/r4/' + name).read().split('\n')
    d = list(difflib.unified_diff(a, b, 'orig', name, n=1, lineterm=''))
    print('\n'.join(d))
    print('--- %d changed hunks' % len([x for x in d if x.startswith('@@')]))
else:
    import capstone
    d, secs, syms = load_elf(mode)
    md = capstone.Cs(capstone.CS_ARCH_ARM64, capstone.CS_MODE_ARM)
    for base, size, nm in syms:
        if not nm.startswith('r4_pass2_8_levels_complex'):
            continue
        text = None
        for s in secs:
            if s['addr'] <= base < s['addr'] + s['size']:
                text = s
                break
        fo = text['offset'] + (base - text['addr'])
        ins = list(md.disasm(d[fo:fo + size], base))
        c = collections.Counter(i.mnemonic for i in ins)
        print('%-38s size=%5d instrs=%5d fp=%4d mem=%4d prfm=%3d fmla=%4d fmul=%4d fadd=%4d fsub=%4d'
              % (nm, size, len(ins),
                 sum(v for k, v in c.items() if k.startswith('f')),
                 sum(v for k, v in c.items() if k.startswith(('ld', 'st'))),
                 c.get('prfm', 0), c.get('fmla', 0), c.get('fmul', 0),
                 c.get('fadd', 0), c.get('fsub', 0)))
