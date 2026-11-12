import json, zlib, sys

# Convert farm results.json.txt lines into submittable JSON by adding the
# checksum object that the mersenne.org manual form requires since 2026.
# Layout mirrors the open mfaktc TF layout, adapted for ECM NF.
LAYOUT = '{exp};ECM;;{curves};{b1};{b2};;;;{prog};{ver};;;{cpu};{ts}'

def add_checksum(line):
    s = LAYOUT.format(exp=line['exponent'], curves=line.get('curves', 1),
                      b1=line['b1'], b2=line['b2'],
                      prog=line['program']['name'], ver=line['program']['version'],
                      cpu='ARM64 Linux', ts=line['timestamp'])
    crc = '%08X' % (zlib.crc32(s.encode()) & 0xFFFFFFFF)
    line['checksum'] = {'version': 1, 'checksum': crc}
    return line

src, dst = sys.argv[1], sys.argv[2]
override_crc = sys.argv[3] if len(sys.argv) > 3 else None
out = []
for l in open(src, encoding='utf-8'):
    l = l.strip()
    if not l:
        continue
    d = add_checksum(json.loads(l))
    if override_crc:
        d['checksum']['checksum'] = override_crc
    out.append(json.dumps(d, separators=(', ', ': ')))
open(dst, 'w', encoding='ascii', newline='\n').write('\n'.join(out) + '\n')
print('wrote', len(out), 'lines to', dst)
