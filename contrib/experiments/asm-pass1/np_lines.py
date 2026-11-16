import sys, struct, collections, bisect, difflib, re

# ---------- ELF object parsing (symbols only) ----------
def load_syms(path):
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

# ---------- .s parsing ----------
DIRECTIVE = re.compile(r'^\s*\.')
def parse_asm(path):
    funcs = {}
    cur = None
    filemap = {}
    fn = None
    fl = 0
    with open(path, 'r', errors='replace') as f:
        for raw in f:
            line = raw.rstrip('\n')
            st = line.strip()
            if not st:
                continue
            if st.startswith('.file'):
                m = re.match(r'\.file\s+(\d+)\s+"(.*)"', st)
                if m:
                    filemap[int(m.group(1))] = m.group(2)
                continue
            if st.startswith('.loc'):
                parts = st.split()
                if len(parts) >= 3:
                    fn = int(parts[1]); fl = int(parts[2])
                continue
            if st.startswith('.type') and 'function' in st:
                name = st.split()[1].split(',')[0].rstrip('@')
                cur = []
                funcs[name] = cur
                continue
            if st.startswith('.size') or st.startswith('.Lfunc_end'):
                if st.startswith('.size'):
                    cur = None
                continue
            if DIRECTIVE.match(st):
                continue
            if cur is None:
                continue
            # instruction: mnemonic + operands
            parts = st.split(None, 1)
            mn = parts[0]
            ops = parts[1] if len(parts) > 1 else ''
            ops = re.split(r'\s*(?://|@)', ops)[0].strip()
            cur.append((mn, ops, filemap.get(fn, str(fn)), fl))
    return funcs

def shape(ops):
    # normalize operands: keep structural shape only
    out = []
    for o in ops.split(','):
        o = o.strip()
        if not o:
            continue
        if o.startswith('#'):
            out.append('#')
        elif o.startswith('['):
            out.append('[')
        elif re.match(r'^[xwsvqdbh]\d+$', o):
            out.append('r')
        elif re.match(r'^[xw]\d+$', o):
            out.append('r')
        else:
            out.append('?')
    return mn_key(ops) + ','.join(out)

def mn_key(ops):
    return ''

def main():
    import capstone
    objpath = sys.argv[1]
    asmpath = sys.argv[2]
    func = sys.argv[3]
    profpath = sys.argv[4] if len(sys.argv) > 4 else None

    d, secs, syms = load_syms(objpath)
    funcs = parse_asm(asmpath)
    if func not in funcs:
        print('func not in .s:', func); return
    sinstrs = funcs[func]
    print('asm instrs for %s: %d' % (func, len(sinstrs)))

    # disassemble object function
    ent = [s for s in syms if s[2] == func]
    if not ent:
        print('func not in obj'); return
    base, size, _ = ent[0]
    text = None
    for s in secs:
        if s['addr'] <= base < s['addr'] + s['size']:
            text = s; break
    fo = text['offset'] + (base - text['addr'])
    code = d[fo:fo + size]
    md = capstone.Cs(capstone.CS_ARCH_ARM64, capstone.CS_MODE_ARM)
    oinstrs = [(i.address, i.mnemonic, i.op_str) for i in md.disasm(code, base)]
    print('obj instrs: %d' % len(oinstrs))

    stok = [mn + ' ' + shape(ops) for mn, ops, f, l in sinstrs]
    otok = [mn + ' ' + shape(ops) for a, mn, ops in oinstrs]
    sm = difflib.SequenceMatcher(None, stok, otok, autojunk=False)
    o2s = {}
    for blk in sm.get_matching_blocks():
        for k in range(blk.size):
            o2s[blk.b + k] = blk.a + k
    print('aligned obj instrs: %d / %d (%.1f%%) ratio=%.3f' % (len(o2s), len(oinstrs), 100.0*len(o2s)/len(oinstrs), sm.ratio()))

    # pc -> obj instr index
    addrs = [a for a, m, o in oinstrs]
    pcs = []
    if profpath:
        raw = open(profpath, 'rb').read()
        t = raw.decode('utf-16') if raw[:2] in (b'\xff\xfe', b'\xfe\xff') else raw.decode('utf-8', 'replace')
        for line in t.splitlines():
            if line.startswith('PC '):
                pcs.append(int(line.split()[1], 16))
    agg = collections.Counter()
    lines = collections.Counter()
    for pc in pcs:
        if not (base <= pc < base + size):
            continue
        i = bisect.bisect_right(addrs, pc) - 1
        agg[i] += 1
        j = o2s.get(i)
        if j is not None:
            mn, ops, f, l = sinstrs[j]
            lines[(f, l)] += 1
    tot = sum(agg.values())
    print('samples: %d' % tot)
    # also: which source files dominate
    fc = collections.Counter()
    for (f, l), c in lines.items():
        fc[f] += c
    print('\nby source file:')
    for f, c in fc.most_common(10):
        print('  %6d %5.1f%%  %s' % (c, 100.0*c/max(tot,1), f))
    # print top source lines with text
    cache = {}
    def getline(f, l):
        key = (f, l)
        if key not in cache:
            try:
                with open(f, 'r', errors='replace') as fh:
                    for i, txt in enumerate(fh, 1):
                        if i == l:
                            cache[key] = txt.rstrip()[:150]
                            break
                    else:
                        cache[key] = ''
            except OSError:
                cache[key] = '<nofile>'
        return cache[key]
    print('\ntop source lines:')
    for (f, l), c in lines.most_common(30):
        print('  %4d %5.1f%%  %s:%d  %s' % (c, 100.0*c/max(tot,1), f.split('\\')[-1].split('/')[-1], l, getline(f, l)))

main()
