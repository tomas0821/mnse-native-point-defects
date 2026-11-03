#!/usr/bin/env python3
"""Decide whether the 149 meV supercell splitting is a band-pairing artefact of folding.

usage: fold_analyse.py uideal_fold.out wzsc_bands.out
 (a) on-path partner only, primitive bands 24-31, index pairing -> should match the ideal-u
     primitive path result (239.0 meV) at this 20 pts/segment sampling
 (b) union of all 18 partners per spin, sorted, index-paired over supercell window 415-558
     -> if ~149 meV, the supercell number is a pairing artefact, not physics
 (c) emulated union eigenvalues vs the actual supercell eigenvalues, bands 1-558
     -> if they agree to ~meV, the pristine supercell IS the folded primitive
"""
import sys
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from splitting import parse

fo, sc = sys.argv[1], sys.argv[2]
d = parse(fo)
(ku, eu), (kd, ed) = d["up"], d["dn"]
assert ku == kd and len(ku) == 41 * 18, len(ku)
NV = 31                                    # valence bands per spin in the primitive
grp = lambda e, j: [e[j * 18 + m] for m in range(18)]

# (a)
a = max(abs(eu[j * 18][n] - ed[j * 18][n]) * 1e3 for j in range(41) for n in range(23, 31))
print(f"(a) on-path partner, bands 24-31, index pairing : {a:.1f} meV")
true_any = max(abs(eu[i][n] - ed[i][n]) * 1e3 for i in range(len(ku)) for n in range(23, 31))
print(f"    max over ALL 738 partner k (true, per-band)   : {true_any:.1f} meV")

# (b)
uni_u = [sorted(x for blk in grp(eu, j) for x in blk[:NV]) for j in range(41)]
uni_d = [sorted(x for blk in grp(ed, j) for x in blk[:NV]) for j in range(41)]
assert all(len(x) == 558 for x in uni_u + uni_d)
b, jb, nb = max((abs(uni_u[j][n] - uni_d[j][n]) * 1e3, j, n)
                for j in range(41) for n in range(414, 558))
print(f"(b) folded union, sorted pairing, window 415-558: {b:.1f} meV  (band {nb+1}, K {jb})")

# (c)
s = parse(sc)
(su, seu), (sd, sed) = s["up"], s["dn"]
assert len(su) == 41
dev = max(abs(seu[j][n] - uni_u[j][n]) for j in range(41) for n in range(558))
dev = max(dev, max(abs(sed[j][n] - uni_d[j][n]) for j in range(41) for n in range(558)))
off = sum(seu[j][n] - uni_u[j][n] for j in range(41) for n in range(558)) / (41 * 558)
print(f"(c) supercell vs emulated eigenvalues, 1-558    : max |dE| = {dev*1e3:.1f} meV, "
      f"mean offset {off*1e3:+.1f} meV")
sc_b = max(abs(seu[j][n] - sed[j][n]) * 1e3 for j in range(41) for n in range(414, 558))
print(f"    actual supercell, window 415-558            : {sc_b:.1f} meV")
