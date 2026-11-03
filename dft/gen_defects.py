#!/usr/bin/env python3
"""
Stage 4 — native point defect supercells.

SUPERCELL CHOICE (settled by dft/build_supercells.py, not by tradition):
the plan's CONVENTIONAL CUBIC sizes are kept. Switching to supercells of the
rhombohedral AFM-II primitive was considered and is NOT better -- efficiency
(min-image distance per atom^(1/3)) is 2.77 for cubic vs 2.82/2.72 for the best
rhombohedral shapes, i.e. a wash, and the cubic cell is easier to reason about
and matches the CrN precedent. The rhombohedral distortion is carried instead as
a uniform +0.289% trigonal strain along [111], applied to the cubic cell.

  rock-salt  2x2x2 conventional = 64 atoms, min image 11.08 A, 3.1% defect conc.
  wurtzite   3x3x2              = 72 atoms, min image 12.81 A, 2.8%  (also the
                                  best shape found in the scan, eff 3.08)

CUTOFF: 120 Ry, matching the chemical-potential references. NOT 70 Ry. The
measured error cancellation (0.02 meV/f.u.) applies to differences between cells
with the SAME atom count. A vacancy changes the atom count, so the cutoff error
on the missing atom (~24 meV/atom at 70 Ry) does not cancel against mu.

Defects, in both polymorphs: V_Se, V_Mn, Mn_Se (Mn on a Se site), Se_Mn.
Neutral only, per the project scope.

    python3 dft/gen_defects.py
"""
import itertools
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "defects")
os.makedirs(OUT, exist_ok=True)

U = 3.7574
A = 5.5378                      # relaxed cubic constant at U=3.76
EPS = -0.00096                  # off-diagonal of the trigonal strain (from the
                                # relaxed AFM-II primitive: F = I + EPS*(J - I))
WZ = np.array([[4.270582592, 0.0, 0.0],
               [-2.135291296, 3.698433014, 0.0],
               [0.0, 0.0, 6.854295216]])
U_WZ = 0.375
PS = {"Mn": "Mn.pbe-spn-kjpaw_psl.0.3.1.UPF", "Se": "Se.pbe-dn-kjpaw_psl.1.0.0.UPF"}
MASS = {"Mn": 54.938044, "Se": 78.971}


def rocksalt_222():
    """conventional 2x2x2, AFM-II, with the trigonal distortion applied."""
    F = np.full((3, 3), EPS) + np.eye(3) * (1 - 2 * EPS - 0.000965 + 0.000965)
    np.fill_diagonal(F, 0.999035)
    cell = (np.eye(3) * 2 * A) @ F
    h = A / 2
    at = []
    for i, j, k in itertools.product(range(4), repeat=3):
        r = np.array([i, j, k]) * h @ F
        if (i + j + k) % 2 == 0:
            at.append(("Mn1" if ((i + j + k) // 2) % 2 == 0 else "Mn2", r))
        else:
            at.append(("Se", r))
    return cell, at


def wurtzite_332():
    a1, a2, a3 = WZ
    frac = [("Mn1", (1/3, 2/3, 0.0)), ("Mn2", (2/3, 1/3, 0.5)),
            ("Se", (1/3, 2/3, U_WZ)), ("Se", (2/3, 1/3, U_WZ + 0.5))]
    at = []
    for i, j, k in itertools.product(range(3), range(3), range(2)):
        t = i * a1 + j * a2 + k * a3
        for lab, f in frac:
            at.append((lab, f[0]*a1 + f[1]*a2 + f[2]*a3 + t))
    cell = np.array([a1 * 3, a2 * 3, a3 * 2])
    return cell, at


def most_central(at, cell, want):
    """pick the site of type `want` closest to the cell centre -- keeps the
    defect as far as possible from the cell faces in the relaxed geometry"""
    c = cell.sum(0) / 2
    idx = [i for i, (l, _) in enumerate(at) if l.startswith(want)]
    return min(idx, key=lambda i: np.linalg.norm(at[i][1] - c))


def make_defect(at, cell, kind):
    at = [list(x) for x in at]
    if kind == "pristine":
        return at, "pristine reference"
    if kind == "V_Se":
        i = most_central(at, cell, "Se"); rm = at.pop(i)
        return at, f"Se vacancy at {np.round(rm[1],3)}"
    if kind == "V_Mn":
        i = most_central(at, cell, "Mn"); rm = at.pop(i)
        return at, f"{rm[0]} vacancy at {np.round(rm[1],3)}"
    if kind == "Mn_Se":                       # Mn sitting on a Se site
        i = most_central(at, cell, "Se"); at[i][0] = "Mn1"
        return at, f"Mn on Se site at {np.round(at[i][1],3)}"
    if kind == "Se_Mn":                       # Se sitting on a Mn site
        i = most_central(at, cell, "Mn"); at[i][0] = "Se"
        return at, f"Se on Mn site at {np.round(at[i][1],3)}"
    raise ValueError(kind)


def write(tag, cell, at, kmesh, note):
    labels = sorted({a[0] for a in at})
    L = ["&CONTROL", "    calculation = 'relax'", f"    prefix      = '{tag}'",
         "    pseudo_dir  = '/home/tomas.rojas_s/pseudo'", f"    outdir      = './tmp_{tag}'",
         "    tprnfor     = .true.", "    tstress     = .true.",
         "    forc_conv_thr = 1.0d-4", "    etot_conv_thr = 1.0d-5", "    nstep = 200", "/",
         "&SYSTEM", "    ibrav       = 0", f"    nat         = {len(at)}",
         f"    ntyp        = {len(labels)}", "    ecutwfc     = 120.0", "    ecutrho     = 960.0",
         "    occupations = 'smearing'", "    smearing    = 'gaussian'", "    degauss     = 0.01",
         "    nspin       = 2"]
    for n, lab in enumerate(labels, start=1):
        m = "+0.5" if lab == "Mn1" else ("-0.5" if lab == "Mn2" else "+0.0")
        L.append(f"    starting_magnetization({n}) = {m}")
    L += ["/", "&ELECTRONS", "    conv_thr        = 1.0d-8", "    mixing_beta     = 0.3",
          "    electron_maxstep = 200", "/", "&IONS", "    ion_dynamics = 'bfgs'", "/",
          "ATOMIC_SPECIES"]
    for lab in labels:
        el = "Mn" if lab.startswith("Mn") else "Se"
        L.append(f"  {lab} {MASS[el]} {PS[el]}")
    if any(l.startswith("Mn") for l in labels):
        L.append("HUBBARD (ortho-atomic)")
        for lab in labels:
            if lab.startswith("Mn"):
                L.append(f"U {lab}-3d {U}")
    L.append("CELL_PARAMETERS angstrom")
    for v in cell:
        L.append("  {:.10f} {:.10f} {:.10f}".format(*v))
    L.append("ATOMIC_POSITIONS angstrom")
    for lab, r in at:
        L.append(f"  {lab} {r[0]:.10f} {r[1]:.10f} {r[2]:.10f}")
    L += ["K_POINTS automatic", "  {} {} {} 0 0 0".format(*kmesh), ""]
    open(os.path.join(OUT, tag + ".in"), "w").write("\n".join(L))
    return len(at), note


KINDS = ["pristine", "V_Se", "V_Mn", "Mn_Se", "Se_Mn"]
print(f"{'tag':22s} {'atoms':>6s}  {'Mn↑/Mn↓/Se':>12s}  note")
made = []
for sysname, builder, km in (("rs", rocksalt_222, (3, 3, 3)),
                             ("wz", wurtzite_332, (3, 3, 3))):
    cell0, at0 = builder()
    for kind in KINDS:
        at, note = make_defect(at0, cell0, kind)
        tag = f"d_{sysname}_{kind}"
        n, _ = write(tag, cell0, at, km, note)
        c = {"Mn1": 0, "Mn2": 0, "Se": 0}
        for l, _ in at:
            c[l] += 1
        print(f"{tag:22s} {n:6d}  {c['Mn1']:4d}/{c['Mn2']:3d}/{c['Se']:3d}  {note}")
        made.append(tag)
print(f"\nwrote {len(made)} inputs to {OUT}")
