#!/usr/bin/env python3
"""Spin splitting |E_up - E_dn| along a k-path from a QE nspin=2 'bands' output.

Pairs bands across the SPIN UP / SPIN DOWN blocks by index and asserts the two
k-lists are identical. Reports peak |split| over a band window, its band and
k-point, the signed peak on each arm (first half / second half of the path),
and the splitting at the path endpoints and midpoint (L1, Gamma, L2 nodes).

usage: splitting.py OUT [--win LO HI] [--all]
"""
import re
import sys

NUM = re.compile(r"-?\d+\.\d+")


def parse(path):
    txt = open(path).read()
    iu = txt.index("------ SPIN UP ------------")
    idn = txt.index("------ SPIN DOWN ----------")
    end = txt.find("the Fermi energy", idn)
    if end < 0:
        end = txt.find("highest occupied", idn)
    if end < 0:
        end = txt.find("Writing", idn)
    blocks = {"up": txt[iu:idn], "dn": txt[idn:end]}
    out = {}
    for s, b in blocks.items():
        parts = re.split(r"\n\s+k =(.*?)\(\s*\d+ PWs\)\s+bands \(ev\):\s*\n", b)
        ks, es = [], []
        for i in range(1, len(parts), 2):
            k = tuple(round(float(x), 4) for x in NUM.findall(parts[i]))
            body = parts[i + 1]
            body = body.split("occupation numbers")[0]
            vals = [float(x) for x in NUM.findall(body)]
            ks.append(k)
            es.append(vals)
        out[s] = (ks, es)
    return out


def main():
    a = sys.argv[1:]
    path = a[0]
    d = parse(path)
    ku, eu = d["up"]
    kd, ed = d["dn"]
    assert ku == kd, "k-lists differ between spin channels"
    nk = len(ku)
    nb = min(len(x) for x in eu + ed)
    assert all(len(x) == len(eu[0]) for x in eu + ed), "ragged band counts"
    lo, hi = 1, nb
    if "--win" in a:
        i = a.index("--win")
        lo, hi = int(a[i + 1]), int(a[i + 2])
    split = [[(eu[k][n] - ed[k][n]) * 1000.0 for n in range(nb)] for k in range(nk)]
    best = max(((abs(split[k][n]), k, n) for k in range(nk) for n in range(lo - 1, hi)))
    mag, kb, nbest = best
    half = nk // 2
    arm1 = [split[k][nbest] for k in range(0, half + 1)]
    arm2 = [split[k][nbest] for k in range(half, nk)]
    a1 = max(arm1, key=abs)
    a2 = max(arm2, key=abs)
    print(f"file      : {path}")
    print(f"nk={nk}  nbnd={nb}  window={lo}-{hi}")
    print(f"peak |split| = {mag:.1f} meV  band {nbest+1}  k = {ku[kb]}  (k index {kb})")
    print(f"arm1 / arm2 signed peak (band {nbest+1}) = {a1:+.1f} / {a2:+.1f} meV")
    nodes = [0, half, nk - 1]
    nv = [max(abs(split[k][n]) for n in range(lo - 1, hi)) for k in nodes]
    print("max |split| in window at L1 / Gamma / L2 = " + " / ".join(f"{x:.2f}" for x in nv) + " meV")
    print(f"   (k = {ku[0]}, {ku[half]}, {ku[nk-1]})")
    if "--all" in a:
        for n in range(lo - 1, hi):
            m = max(abs(split[k][n]) for k in range(nk))
            print(f"  band {n+1:4d}  max|split| {m:8.1f} meV")


if __name__ == "__main__":
    main()
