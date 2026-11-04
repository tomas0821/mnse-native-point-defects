#!/usr/bin/env python3
"""Kohn-Sham levels of the defect cells relative to the host valence-band maximum, and smearing term.

From the final eigenvalue block of each 120 Ry output: E_F, highest level below and lowest level above it
per spin (these outputs do not print occupations; constrained cells have one E_F per spin), aligned between cell and host through the k-weighted Se-4s
centroid, and -TS. Levels are Kohn-Sham energies on the 3x3x3 grid, not transition levels. In gapped cells
E_F only reflects where Gaussian neutrality places it between HOMO and LUMO.

usage: ks_levels.py   (prints the table; import rows() for the SI table)
"""
import os
import re

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..")
RY = 13.605693123
CELLS = {
    "rs": {"host": "prod/d_rs_pristine_s120.out", "V_Mn": "prod/d_rs_V_Mn_s120.out",
           "V_Se": "prod/d_rs_V_Se_s120.out", "Se_Mn": "prod/d_rs_Se_Mn_s120.out",
           "Mn_Se": "anti/rsMnSe_M5_s120.out"},
    "wz": {"host": "prod/d_wz_pristine_s120.out", "V_Mn": "prod/d_wz_V_Mn_s120.out",
           "V_Se": "prod/d_wz_V_Se_s120.out", "Se_Mn": "prod/d_wz_Se_Mn_s120.out",
           "Mn_Se": "prod/d_wz_Mn_Se_s120.out"},
}
WIN = {"rs": (-6.0, 1.0), "wz": (-9.0, -3.0)}   # Se-4s manifold of each polymorph (eV)


def parse(path):
    t = open(path, errors="ignore").read()
    wk = [float(x) for x in re.findall(r"k\(\s*\d+\) = \(.*?\), wk =\s+([\d.]+)", t)]
    idx = t.rfind("End of self-consistent calculation")
    blk = t[idx:]
    up, dn = blk.index("SPIN UP"), blk.index("SPIN DOWN")
    end = min(x for x in (blk.find("the Fermi energy", dn), blk.find("the spin up/dw Fermi", dn),
                          blk.find("highest occupied", dn)) if x > 0)
    out = {}
    for s, b in (("up", blk[up:dn]), ("dn", blk[dn:end])):
        E = []
        for m in re.finditer(r"bands \(ev\):\s*\n(.*?)(?=\n\s*k =|\Z)", b, re.S):
            E.append([float(x) for x in re.findall(r"-?\d+\.\d+", m.group(1))])
        out[s] = np.array(E)
    nk = out["up"].shape[0]
    ef = re.search(r"the Fermi energy is\s+(-?[\d.]+)", blk)
    ef2 = re.search(r"the spin up/dw Fermi energies are\s+(-?[\d.]+)\s+(-?[\d.]+)", blk)
    hl = re.search(r"highest occupied, lowest unoccupied level \(ev\):\s+(-?[\d.]+)\s+(-?[\d.]+)", blk)
    if ef:
        f = (float(ef.group(1)),) * 2
    elif ef2:
        f = (float(ef2.group(1)), float(ef2.group(2)))
    else:
        mid = (float(hl.group(1)) + float(hl.group(2))) / 2
        f = (mid, mid)
    ts = re.findall(r"smearing contrib. \(-TS\)\s+=\s+(-?[\d.]+) Ry", t)
    mag = re.findall(r"total magnetization\s+=\s+(-?[\d.]+)", t)
    return dict(E=out, w=np.array(wk[:nk]), ef=f, single=bool(ef),
                ts=float(ts[-1]) * RY * 1e3 if ts else 0.0, M=float(mag[-1]) if mag else 0.0)


def se4s(r, win):
    a = b = 0.0
    for s in ("up", "dn"):
        E = r["E"][s]
        w = np.repeat(r["w"][:, None], E.shape[1], axis=1)
        m = (E > win[0]) & (E < win[1])
        a += (E[m] * w[m]).sum()
        b += w[m].sum()
    return a / b


def rows():
    res = {}
    for p, cells in CELLS.items():
        h = parse(os.path.join(D, cells["host"]))
        c0 = se4s(h, WIN[p])
        vbm = max(h["E"][s][h["E"][s] < h["ef"][0]].max() for s in ("up", "dn"))
        for d, f in cells.items():
            r = parse(os.path.join(D, f))
            sh = se4s(r, WIN[p]) - c0
            lev = {}
            for i, s in enumerate(("up", "dn")):
                E = r["E"][s]
                lev[s] = (E[E < r["ef"][i]].max() - sh - vbm, E[E >= r["ef"][i]].min() - sh - vbm)
            ef = r["ef"][0] - sh - vbm if r["single"] else None
            res[(p, d)] = dict(ef=ef, homo=(lev["up"][0], lev["dn"][0]), lumo=(lev["up"][1], lev["dn"][1]),
                               ts=r["ts"], M=r["M"])
    return res


if __name__ == "__main__":
    for k, v in rows().items():
        ef = "  --  " if v["ef"] is None else f"{v['ef']:+.2f}"
        print(f"{k[0]} {k[1]:6s} E_F {ef} HOMO {v['homo'][0]:+.2f}/{v['homo'][1]:+.2f} "
              f"LUMO {v['lumo'][0]:+.2f}/{v['lumo'][1]:+.2f} -TS {v['ts']:7.1f} meV M {v['M']:+.2f}")
