#!/usr/bin/env python3
"""
QE relax/scf output  ->  publication figures.

Writes, for each structure:
  <tag>.cif    clean structure, opens in anything
  <tag>.vesta  VESTA scene with SPIN ARROWS already attached as vector properties,
               plus per-species colours -- no manual vector entry in the GUI

Why this route: VESTA is the field standard for solid-state structure figures and
supports per-atom vectors natively, but it is GUI-driven and there are ten defect
cells. This script does the mechanical part so VESTA is only used for the final
aesthetic pass (camera, lighting, polyhedra, export).

ASE's own espresso-out reader asserts on spin-polarised output
(`AssertionError: ((2,0), 6)` -- it expects one eigenvalue set per k-point), so the
geometry and moments are parsed directly here and handed to ASE clean.

    python3 dft/qe2vesta.py <qe_output> [more outputs...] -o figures/structures
"""
import os
import re
import sys

import numpy as np

# element colours (VESTA-ish) and radii; Mn split by spin sign
COL = {"Mn_up": (0.24, 0.37, 0.64), "Mn_dn": (0.76, 0.27, 0.18), "Se": (0.62, 0.65, 0.71)}
RAD = {"Mn_up": 0.60, "Mn_dn": 0.60, "Se": 0.45}


def parse_qe(path):
    """last geometry + per-site moments from a pw.x relax/scf output"""
    txt = open(path, errors="ignore").read()

    # --- cell: prefer the final relaxed one, else the input cell
    cells = re.findall(r"CELL_PARAMETERS \((?:angstrom|bohr)\)\n(.*?)\n\s*\n", txt, re.S)
    if cells:
        cell = np.array([[float(x) for x in r.split()] for r in cells[-1].strip().splitlines()])
        if "bohr" in txt.split("CELL_PARAMETERS")[1][:20]:
            cell *= 0.529177
    else:
        m = re.search(r"crystal axes: \(cart\. coord\. in units of alat\)\n(.*?)\n\s*\n", txt, re.S)
        alat = float(re.search(r"lattice parameter \(alat\)\s*=\s*([\d.]+)", txt).group(1)) * 0.529177
        cell = np.array([[float(x) for x in r.split()[3:6]] for r in m.group(1).splitlines()]) * alat

    # --- positions: last ATOMIC_POSITIONS block
    blocks = re.findall(r"ATOMIC_POSITIONS \((\w+)\)\n(.*?)(?=\n\s*\n|\nEnd final)", txt, re.S)
    if not blocks:
        # a plain scf never reprints ATOMIC_POSITIONS -- fall back to the header
        # table, which is in alat units
        m = re.search(r"site n\.\s+atom\s+positions \(alat units\)\n(.*?)\n\s*\n", txt, re.S)
        if not m:
            raise ValueError(f"no atomic positions found in {path}")
        alat = float(re.search(r"lattice parameter \(alat\)\s*=\s*([\d.]+)", txt).group(1)) * 0.529177
        labels, pos = [], []
        for ln in m.group(1).splitlines():
            g = re.match(r"\s*\d+\s+(\S+)\s+tau\(\s*\d+\)\s*=\s*\(\s*(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)", ln)
            if g:
                labels.append(g.group(1))
                pos.append([float(g.group(i)) * alat for i in (2, 3, 4)])
        pos = np.array(pos)
        mag = np.zeros(len(labels))
        mb = re.findall(r"Magnetic moment per site.*?\n((?:\s*atom\s+\d+.*\n)+)", txt)
        if mb:
            vals = re.findall(r"magn=\s*(-?\d+\.\d+)", mb[-1])
            if len(vals) == len(labels):
                mag = np.array([float(v) for v in vals])
        return cell, labels, pos, mag
    unit, body = blocks[-1]
    labels, pos = [], []
    for ln in body.strip().splitlines():
        p = ln.split()
        if len(p) < 4:
            continue
        labels.append(p[0])
        pos.append([float(x) for x in p[1:4]])
    pos = np.array(pos)
    if unit == "crystal":
        pos = pos @ cell
    elif unit == "bohr":
        pos *= 0.529177

    # --- per-site magnetic moments (last printed block)
    mag = np.zeros(len(labels))
    mb = re.findall(r"Magnetic moment per site.*?\n((?:\s*atom\s+\d+.*\n)+)", txt)
    if mb:
        vals = re.findall(r"magn=\s*(-?\d+\.\d+)", mb[-1])
        if len(vals) == len(labels):
            mag = np.array([float(v) for v in vals])
    return cell, labels, pos, mag


def species(label, m):
    if label.startswith("Mn"):
        return "Mn_up" if m >= 0 else "Mn_dn"
    return "Se"


def write_vesta(path, cell, labels, pos, mag, title):
    """minimal but valid .vesta scene, with spin vectors in VECTR/VECTT"""
    import math
    a, b, c = [np.linalg.norm(v) for v in cell]
    al = math.degrees(math.acos(cell[1] @ cell[2] / (b * c)))
    be = math.degrees(math.acos(cell[0] @ cell[2] / (a * c)))
    ga = math.degrees(math.acos(cell[0] @ cell[1] / (a * b)))
    inv = np.linalg.inv(cell)
    frac = pos @ inv

    L = ["#VESTA_FORMAT_VERSION 3.5.0", "", "CRYSTAL", "", "TITLE", title, "",
         "GROUP", "1 1 P 1", "SYMOP",
         " 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1   1",
         " 0.000000  0.000000  0.000000  0  0  0   0  0  0   0  0  0   0",
         "TRANM 0",
         " 0.000000  0.000000  0.000000  1  0  0   0  1  0   0  0  1",
         "LTRANSL", " -1", " 0.000000  0.000000  0.000000  0.000000  0.000000  0.000000",
         "LORIENT", " -1   0   0   0   0",
         " 1.000000  0.000000  0.000000  1.000000  0.000000  0.000000",
         " 0.000000  0.000000  1.000000  0.000000  0.000000  1.000000",
         "LMATRIX",
         " 1.000000  0.000000  0.000000  0.000000", " 0.000000  1.000000  0.000000  0.000000",
         " 0.000000  0.000000  1.000000  0.000000", " 0.000000  0.000000  0.000000  1.000000",
         " 0.000000  0.000000  0.000000",
         "CELLP", f" {a:.6f} {b:.6f} {c:.6f} {al:.6f} {be:.6f} {ga:.6f}",
         " 0.000000 0.000000 0.000000 0.000000 0.000000 0.000000", "STRUC"]
    for i, (lab, f, m) in enumerate(zip(labels, frac, mag), start=1):
        el = "Mn" if lab.startswith("Mn") else "Se"
        sp = species(lab, m)
        L.append(f" {i} {el} {sp}{i} 1.0000 {f[0]:.6f} {f[1]:.6f} {f[2]:.6f} 1a 1")
        L.append("                            0.000000   0.000000   0.000000  0.00")
    L += [" 0 0 0 0 0 0 0", "THERI 0"]
    for i, lab in enumerate(labels, start=1):
        L.append(f" {i} {lab}{i} 1.000000")
    L += [" 0 0 0", "SHAPE", "  0       0       0       0   0.000000  0   192   192   192   192",
          "BOUND", "       0        1         0        1         0        1",
          " 0   0   0   0  0", "SBOND", " 0 0 0 0", "SITET"]
    for i, (lab, m) in enumerate(zip(labels, mag), start=1):
        sp = species(lab, m)
        r, g, bl = COL[sp]
        L.append(f" {i} {sp}{i} {RAD[sp]:.4f} {int(r*255)} {int(g*255)} {int(bl*255)} "
                 f"{int(r*255)} {int(g*255)} {int(bl*255)} 204 0")
    L += [" 0 0 0 0 0 0", "VECTR"]
    nv = 0
    for i, m in enumerate(mag, start=1):
        if abs(m) < 0.2:            # only magnetic sites get an arrow
            continue
        nv += 1
        L.append(f" {nv}    0.000000    0.000000 {m*0.35:.6f} 0")   # arrow along +/- c
        L.append(f"    {i}   0    0    0    0")
        L.append(" 0 0 0 0 0")
    L += [" 0 0 0 0 0", "VECTT"]
    for k in range(1, nv + 1):
        L.append(f" {k}  0.500 30 30 30 1")     # radius, RGB (dark), penetrate flag
    L += [" 0 0 0 0 0", "SPLAN", "  0   0   0   0", "LBLAT", " -1", "LBLSP", " -1",
          "DLATM", " -1", "DLBND", " -1", "DLPLY", " -1",
          "ATOMT"]
    for n, sp in enumerate(("Mn_up", "Mn_dn", "Se"), start=1):
        r, g, bl = COL[sp]
        L.append(f"  {n} {sp} {RAD[sp]:.4f} {int(r*255)} {int(g*255)} {int(bl*255)} "
                 f"{int(r*255)} {int(g*255)} {int(bl*255)} 204")
    L += ["  0 0 0 0 0 0", "SCENE",
          "-0.850000  0.000000  0.000000  0.526650",
          " 0.000000  0.000000  0.000000", " 0.000000", " 1.000000",
          "HBOND 0 2", "", "STYLE", "DISPF 37753794", "MODEL   0  1  0",
          "SURFS   0  1  1", "SECTP", "  1  5.00000E+00  5.00000E+00  0.00000E+00",
          "FORMS   0  1", "ATOMS   0  0  1", "BONDS   1", "POLYS   1",
          "TEXCL   1", "VECTS 1.000000", "FORMP", "  1  1.0   0   0   0",
          "ATOMP", "  24  24  0  50  2.0  0", "BONDP",
          "  1  16  0.250  2.000 127 127 127", "POLYP",
          " 204 1  1.000 180 180 180", ""]
    open(path, "w").write("\n".join(L))
    return nv


def main(argv):
    out = "figures/structures"
    if "-o" in argv:
        i = argv.index("-o"); out = argv[i + 1]; argv = argv[:i] + argv[i + 2:]
    os.makedirs(out, exist_ok=True)
    from ase import Atoms
    from ase.io import write as ase_write
    for src in argv:
        tag = os.path.basename(src).replace(".out", "")
        cell, labels, pos, mag = parse_qe(src)
        syms = ["Mn" if l.startswith("Mn") else "Se" for l in labels]
        at = Atoms(symbols=syms, positions=pos, cell=cell, pbc=True)
        at.set_initial_magnetic_moments(mag)
        ase_write(os.path.join(out, tag + ".cif"), at)
        nv = write_vesta(os.path.join(out, tag + ".vesta"), cell, labels, pos, mag, tag)
        up = sum(1 for m in mag if m > 0.2); dn = sum(1 for m in mag if m < -0.2)
        print(f"  {tag:22s} {len(at):3d} atoms  Mn↑{up:3d} Mn↓{dn:3d}  "
              f"{nv} spin arrows  |m|max={np.abs(mag).max():.2f}")


if __name__ == "__main__":
    main(sys.argv[1:])
