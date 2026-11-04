#!/usr/bin/env python3
"""Wurtzite MnSe magnetic orders in the doubled cell 2a1 x a2 x c (8 atoms, 4 Mn).

Tests the stripe orders of Mavani et al., Phys. Rev. Materials 10, 064408 (2026), against the
A-type order used throughout this project, at this project's settings (U = 3.7574 eV, Dudarev,
ortho-atomic, 120/960 Ry). A-type and FM are run in the same doubled cell as controls: A-type must
reproduce the 4-atom-cell energy, which measures the k-grid and cell error of the comparison.

Each layer holds 2 Mn. Spins are assigned per site; the script counts, for every Mn tetrahedron of
the hcp Mn sublattice, the number of up and down spins:
  AFMs-alpha: every tetrahedron 2:2   AFMs-beta: tetrahedra 3:1   A-type: tetrahedra 3:1

Usage: gen_stripe.py   (writes the inputs and prints the tetrahedron census; no calculation here)
"""
import itertools
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# Relaxed 4-atom cell at U = 3.7574 eV (usweep/wz_U37574.out, vc-relax final coordinates, A-type).
A = np.array([[4.270583028, 0.0, 0.0],
              [-2.135291514, 3.698433391, 0.0],
              [0.0, 0.0, 6.854295104]])
MN = np.array([[0.0, 2.4656222608, -0.0096777392],
               [2.1352915140, 1.2328111304, 3.4174697914]])
SE = np.array([[0.0, 2.4656222608, 2.5800391630],
               [2.1352915140, 1.2328111304, 6.0071852168]])

CELL = A * np.array([[2], [1], [1]])
# Mn sites in the doubled cell: (layer, chain)
SITES = [MN[0], MN[0] + A[0], MN[1], MN[1] + A[0]]
SE_SITES = [SE[0], SE[0] + A[0], SE[1], SE[1] + A[0]]

# Spin of each of the four Mn sites: layer A chains 0,1 then layer B chains 0,1.
ORDERS = {
    "Atype": (+1, +1, -1, -1),
    "FM": (+1, +1, +1, +1),
    "AFMsb": (+1, -1, +1, -1),  # Mavani AFMs-beta: tetrahedra 3:1
    "AFMsa": (+1, -1, -1, +1),  # Mavani AFMs-alpha: tetrahedra 2:2
}

INPUT = """&CONTROL
    calculation = '{calc}'
    prefix      = '{name}'
    pseudo_dir  = '/home/tomas.rojas_s/pseudo'
    outdir      = './tmp_{name}'
    tprnfor     = .true.
    tstress     = {tstress}
    forc_conv_thr = 1.0d-4
    etot_conv_thr = 1.0d-5
    nstep       = 100
/
&SYSTEM
    ibrav       = 0
    nat         = 8
    ntyp        = {ntyp}
    ecutwfc     = 120.0
    ecutrho     = 960.0
    occupations = 'smearing'
    smearing    = 'gaussian'
    degauss     = 0.01
    nspin       = 2
{starts}
/
&ELECTRONS
    conv_thr        = 1.0d-8
    mixing_beta     = 0.3
    electron_maxstep = 200
/
{extra}ATOMIC_SPECIES
{species}
HUBBARD (ortho-atomic)
{hubbard}
CELL_PARAMETERS angstrom
{cell}
ATOMIC_POSITIONS angstrom
{pos}
K_POINTS automatic
  5 10 6 0 0 0
"""

RELAX_EXTRA = """&IONS
    ion_dynamics = 'bfgs'
/
&CELL
    cell_dynamics = 'bfgs'
    press_conv_thr = 0.5
/
"""


def tetrahedra():
    """Mn tetrahedra of the hcp sublattice: every set of 4 mutually nearest-neighbour Mn sites,
    found on a large periodic image set and reduced to doubled-cell site indices."""
    images = []
    for n in itertools.product(range(-2, 3), range(-2, 3), range(-1, 2)):
        for i, s in enumerate(SITES):
            images.append((i, s + np.dot(n, CELL)))
    pos = np.array([p for _, p in images])
    idx = [i for i, _ in images]
    home = [k for k, (i, p) in enumerate(images)
            if np.linalg.norm(p - SITES[i]) < 1e-6]
    d = np.linalg.norm(pos[:, None] - pos[None, :], axis=2)
    cut = 4.45  # between the nn (4.27-4.28 A) and next shell (6.04 A)
    found = {}
    for h in home:
        nb = [k for k in range(len(pos)) if 0 < d[h, k] < cut]
        for a, b, c in itertools.combinations(nb, 3):
            if d[a, b] < cut and d[a, c] < cut and d[b, c] < cut:
                key = tuple(sorted(np.round(pos[[h, a, b, c]].mean(0) % 1e3, 3).tolist()))
                found.setdefault(key, sorted((idx[h], idx[a], idx[b], idx[c])))
    return list(found.values())


def census(order, tets):
    out = {}
    for t in tets:
        up = sum(order[i] > 0 for i in t)
        key = f"{max(up, 4 - up)}:{min(up, 4 - up)}"
        out[key] = out.get(key, 0) + 1
    return out


def write(name, order, relax):
    mags = sorted(set(order), reverse=True)
    lab = {+1: "Mn1", -1: "Mn2"}
    used = [lab[m] for m in mags]
    species = [f"  {l} 54.938044 Mn.pbe-spn-kjpaw_psl.0.3.1.UPF" for l in used]
    species.append("  Se 78.971 Se.pbe-dn-kjpaw_psl.1.0.0.UPF")
    starts = [f"    starting_magnetization({i + 1}) = {0.5 * m:+.1f}" for i, m in enumerate(mags)]
    starts.append(f"    starting_magnetization({len(mags) + 1}) = +0.0")
    pos = [f"  {lab[o]} {s[0]:.10f} {s[1]:.10f} {s[2]:.10f}" for o, s in zip(order, SITES)]
    pos += [f"  Se {s[0]:.10f} {s[1]:.10f} {s[2]:.10f}" for s in SE_SITES]
    text = INPUT.format(
        calc="vc-relax" if relax else "scf",
        name=name,
        tstress=".true." if relax else ".false.",
        ntyp=len(used) + 1,
        starts="\n".join(starts),
        extra=RELAX_EXTRA if relax else "",
        species="\n".join(species),
        hubbard="\n".join(f"U {l}-3d 3.7574" for l in used),
        cell="\n".join("  " + " ".join(f"{x:.10f}" for x in v) for v in CELL),
        pos="\n".join(pos),
    )
    with open(os.path.join(HERE, name + ".in"), "w") as f:
        f.write(text)


def main():
    tets = tetrahedra()
    print(f"{len(tets)} distinct Mn tetrahedra in the doubled cell")
    for name, order in ORDERS.items():
        print(f"  {name:8s} spins {order}: tetrahedra {census(order, tets)}")
        write(name, order, relax=False)
        if name in ("Atype", "AFMsa"):
            write(name + "_vc", order, relax=True)


if __name__ == "__main__":
    main()
