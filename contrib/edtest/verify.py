#!/usr/bin/env python3
"""Verify edtest output against independent big-int math.

edtest prints one value per line:  TAG hexvalue
Recompute every operation mod N = 2^n - 1 and report PASS/FAIL per test.
"""
import sys

def load(path):
    vals = {}
    order = []
    cur_n = None
    for ln in open(path, encoding='ascii'):
        ln = ln.strip()
        if ln.startswith('==='):
            cur_n = int(ln.split()[1].split('=')[1])
            vals[cur_n] = {}
            order.append(cur_n)
            continue
        if not ln or cur_n is None:
            continue
        tag, hexs = ln.split(None, 1)
        vals[cur_n][tag] = int(hexs, 16)
    return order, vals

def report(name, got, want):
    ok = (got == want % (2**0)) and (got == want)  # want already reduced
    if got == want:
        print("  PASS %s" % name)
        return True
    # show first differing bit for diagnosis
    x = got ^ want
    bit = x.bit_length() - 1
    print("  FAIL %s: first diff bit %d (got %d bits, want %d bits)" %
          (name, bit, got.bit_length(), want.bit_length()))
    return False

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'edtest.out'
    order, vals = load(path)
    allok = True
    for n in order:
        N = (1 << n) - 1
        v = vals[n]
        print("n=%d: %d values" % (n, len(v)))
        ok = True
        if 'BASIC_a' in v:
            ok &= report('basic mul', v['BASIC_ab'], v['BASIC_a'] * v['BASIC_b'] % N)
        if 'MMA5_a' in v:
            s1, s2, s3, s4 = v['MMA5_a'], v['MMA5_b'], v['MMA5_c'], v['MMA5_d']
            want = (s1*s2 + s3*s4) % N
            ok &= report('mulmuladd5', v['MMA5_abcd'], want)
            ok &= report('mulmuladd5+SNF', v['MMA5S_abcd'], want)
            ok &= report('submul4', v['SM4_amb'], (s1-s2)*s3 % N)
            ok &= report('addsub4o add', v['AS4_add'], (s1+s2) % N)
            ok &= report('addsub4o sub', v['AS4_sub'], (s1-s2) % N)
        if 'DBL_XYZ_x' in v:
            X, Y, Z = v['DBL_XYZ_x'], v['DBL_XYZ_y'], v['DBL_XYZ_z']
            C = 2 * Z * Z % N
            B = 2 * X * X % N
            E = 2 * X * Y % N
            G = (X * X + Y * Y) % N
            F = (C - G) % N
            H = (G - B) % N
            ok &= report('dbl_xyz X3', v['DBL_XYZ_X3'], F * E % N)
            ok &= report('dbl_xyz Y3', v['DBL_XYZ_Y3'], H * G % N)
            ok &= report('dbl_xyz Z3', v['DBL_XYZ_Z3'], F * G % N)
        if 'DBLF_x' in v:
            X, Y, Z = v['DBLF_x'], v['DBLF_y'], v['DBLF_z']
            C = 2 * Z * Z % N
            B = 2 * X * X % N
            E = 2 * X * Y % N
            G = (X * X + Y * Y) % N
            F = (C - G) % N
            H = (G - B) % N
            ok &= report('dblf X3', v['DBLF_X3'], F * E % N)
            ok &= report('dblf Y3', v['DBLF_Y3'], H * G % N)
            ok &= report('dblf Z3', v['DBLF_Z3'], F * G % N)
            ok &= report('dblf T3', v['DBLF_T3'], H * E % N)
        for pfx in ('ADDX', 'ADDS'):
            x1 = pfx + '_x1'
            if x1 not in v:
                continue
            X1, Y1, T1, Z1 = v[pfx+'_x1'], v[pfx+'_y1'], v[pfx+'_t1'], v[pfx+'_z1']
            X2, Y2, T2, Z2 = v[pfx+'_x2'], v[pfx+'_y2'], v[pfx+'_t2'], v[pfx+'_z2']
            D = T1 * Z2 % N
            C = Z1 * T2 % N
            E = (D + C) % N
            H = (D - C) % N
            F = (X1 * Y2 - Y1 * X2) % N
            G = (Y1 * Y2 + X1 * X2) % N
            ok &= report(pfx+' X3', v[pfx+'_X3'], E * F % N)
            ok &= report(pfx+' Y3', v[pfx+'_Y3'], G * H % N)
            ok &= report(pfx+' Z3', v[pfx+'_Z3'], F * G % N)
        allok &= ok
        print("  => n=%d %s" % (n, 'ALL PASS' if ok else 'FAILURES'))
    print('OVERALL:', 'PASS' if allok else 'FAIL')

main()
