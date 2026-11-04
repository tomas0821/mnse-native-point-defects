#!/usr/bin/env python3
"""Does the magnetic order of the wurtzite host change the defect formation energies?

Single points at fixed geometry: the 4 x 2 x 2 supercell of the primitive cell (64 atoms, relaxed
A-type lattice and u = 0.3778, no atomic relaxation), which carries both the A-type order and the
stripe order AFMs-alpha of Mavani et al. The 3 x 3 x 2 production cell cannot carry the stripe order,
because it has an odd repeat along a1. Both orders use the same cell, geometry, grid and cutoff, so the
difference isolates the magnetic contribution to E_f:

    dE_f(D) = [E(D) - E(host) - n_Mn * E(host)/32]_stripe - [same]_Atype

Usage: gen_stripe_defects.py  (writes inputs; prints the Mn neighbour census of each defect site)
"""
import itertools
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

A = np.array([[4.270583028, 0.0, 0.0],
              [-2.135291514, 3.698433391, 0.0],
              [0.0, 0.0, 6.854295104]])
MN = np.array([[0.0, 2.4656222608, -0.0096777392],
               [2.1352915140, 1.2328111304, 3.4174697914]])
SE = np.array([[0.0, 2.4656222608, 2.5800391630],
               [2.1352915140, 1.2328111304, 6.0071852168]])
REP = (4, 2, 2)
CELL = A * np.array(REP)[:, None]


def spin(order, i, layer):
    """Spin of the Mn at repeat i along a1 in layer 0 (A) or 1 (B)."""
    sub = +1 if layer == 0 else -1
    return sub if order == "Atype" else sub * (-1) ** i


def build(order):
    atoms = []
    for i, j, k in itertools.product(*map(range, REP)):
        t = i * A[0] + j * A[1] + k * A[2]
        for layer in (0, 1):
            atoms.append({"el": "Mn", "pos": MN[layer] + t, "s": spin(order, i, layer),
                          "key": ("Mn", i, j, k, layer)})
    for i, j, k in itertools.product(*map(range, REP)):
        t = i * A[0] + j * A[1] + k * A[2]
        for layer in (0, 1):
            atoms.append({"el": "Se", "pos": SE[layer] + t, "s": 0, "key": ("Se", i, j, k, layer)})
    return atoms


def neighbours(atoms, center, el, cut=2.9):
    inv = np.linalg.inv(CELL)
    out = []
    for a in atoms:
        if a["el"] != el:
            continue
        d = a["pos"] - center
        d -= np.round(d @ inv) @ CELL
        if 0.1 < np.linalg.norm(d) < cut:
            out.append(a)
    return out


TEMPLATE = """&CONTROL
    calculation = 'scf'
    prefix      = '{name}'
    pseudo_dir  = '/home/tomas.rojas_s/pseudo'
    outdir      = './tmp_{name}'
    tprnfor     = .false.
    tstress     = .false.
    disk_io     = 'low'
/
&SYSTEM
    ibrav       = 0
    nat         = {nat}
    ntyp        = 3
    ecutwfc     = 70.0
    ecutrho     = 560.0
    occupations = 'smearing'
    smearing    = 'gaussian'
    degauss     = 0.01
    nspin       = 2
    starting_magnetization(1) = +0.5
    starting_magnetization(2) = -0.5
    starting_magnetization(3) = +0.0
/
&ELECTRONS
    conv_thr        = 1.0d-8
    mixing_beta     = 0.2
    electron_maxstep = 300
/
ATOMIC_SPECIES
  Mn1 54.938044 Mn.pbe-spn-kjpaw_psl.0.3.1.UPF
  Mn2 54.938044 Mn.pbe-spn-kjpaw_psl.0.3.1.UPF
  Se 78.971 Se.pbe-dn-kjpaw_psl.1.0.0.UPF
HUBBARD (ortho-atomic)
U Mn1-3d 3.7574
U Mn2-3d 3.7574
CELL_PARAMETERS angstrom
{cell}
ATOMIC_POSITIONS angstrom
{pos}
K_POINTS automatic
  2 4 3 0 0 0
"""


def write(name, atoms):
    lab = {+1: "Mn1", -1: "Mn2", 0: "Se"}
    pos = "\n".join(f"  {lab[a['s']]} {a['pos'][0]:.10f} {a['pos'][1]:.10f} {a['pos'][2]:.10f}"
                    for a in atoms)
    cell = "\n".join("  " + " ".join(f"{x:.10f}" for x in v) for v in CELL)
    with open(os.path.join(HERE, name + ".in"), "w") as f:
        f.write(TEMPLATE.format(name=name, nat=len(atoms), cell=cell, pos=pos))


def main():
    for order in ("Atype", "AFMsa"):
        atoms = build(order)
        n_up = sum(a["s"] > 0 for a in atoms)
        print(f"{order}: {len(atoms)} atoms, Mn up {n_up}, down {sum(a['s'] < 0 for a in atoms)}")
        write(f"sd_{order}_host", atoms)

        vmn = ("Mn", 0, 0, 0, 0)
        vse = ("Se", 0, 0, 0, 0)
        site = {a["key"]: a for a in atoms}
        for nm, key in (("Mn site", vmn), ("Se site", vse)):
            nb = neighbours(atoms, site[key]["pos"], "Mn", cut=2.9 if key[0] == "Se" else 4.45)
            print(f"  {nm}: {len(nb)} Mn neighbours, up {sum(a['s'] > 0 for a in nb)}, "
                  f"down {sum(a['s'] < 0 for a in nb)}")

        write(f"sd_{order}_VMn", [a for a in atoms if a["key"] != vmn])
        write(f"sd_{order}_VSe", [a for a in atoms if a["key"] != vse])
        sm = []
        for a in atoms:
            if a["key"] == vmn:
                sm.append(dict(a, el="Se", s=0))
            else:
                sm.append(a)
        write(f"sd_{order}_SeMn", sm)


if __name__ == "__main__":
    main()
