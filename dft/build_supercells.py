#!/usr/bin/env python3
"""
Stage 3 — defect supercell design, and the geometry check that should precede it.

The original plan specified rock-salt 2x2x2 (64 atoms) and 3x3x3 (216) -- but
those are CONVENTIONAL CUBIC cells, written before we knew the AFM-II ground
state is rhombohedrally distorted (+0.129 deg). Two candidate constructions:

  (A) conventional cubic supercells  -- the plan's sizes, but the distortion has
      to be imposed on a cubic box by hand
  (B) supercells of the RELAXED RHOMBOHEDRAL AFM-II primitive cell (4 atoms) --
      distortion and type-II order come along automatically, by construction

What actually controls the spurious defect-defect interaction is not the atom
count but the MINIMUM IMAGE DISTANCE: how far a defect sits from its own
periodic images. This script computes that for both families so the size choice
is made on the number that matters rather than on tradition.

    python3 dft/build_supercells.py
"""
import itertools
import math

import numpy as np

# ---- relaxed cells at U = 3.76 (from s3_rs_AFMII_prim and b_wz_relax) --------
RS_PRIM = np.array([[-2.768886205,  2.768886205,  0.000000000],
                    [-0.000000000, -2.768886205,  2.768886205],
                    [ 2.758257176,  2.758257176,  5.527143381]])
A_RS_CUBIC = 5.5378          # implied conventional cubic constant
WZ = np.array([[ 4.270582592,  0.0,          0.0        ],
               [-2.135291296,  3.698433014,  0.0        ],
               [ 0.0,          0.0,          6.854295216]])


def min_image(cell, n):
    """Shortest lattice vector of the n1 x n2 x n3 supercell = the closest a
    defect ever comes to one of its own periodic images."""
    S = np.array([cell[i] * n[i] for i in range(3)])
    best = 1e9
    for c in itertools.product(range(-2, 3), repeat=3):
        if c == (0, 0, 0):
            continue
        best = min(best, np.linalg.norm(np.array(c) @ S))
    return best


def report(name, cell, atoms_per_cell, mn_per_cell, sizes):
    print(f"\n{name}")
    print(f"  {'supercell':>10s} {'atoms':>7s} {'Mn sites':>9s} {'min image (A)':>14s} "
          f"{'defect conc.':>13s}")
    for n in sizes:
        mult = n[0] * n[1] * n[2]
        nat = atoms_per_cell * mult
        nmn = mn_per_cell * mult
        d = min_image(cell, n)
        print(f"  {str(n):>10s} {nat:7d} {nmn:9d} {d:14.3f} {1/nmn:12.3%}")


print("=" * 74)
print("DEFECT-DEFECT SEPARATION: what actually limits supercell size")
print("=" * 74)

# (B) rhombohedral AFM-II primitive: 4 atoms (2 Mn + 2 Se), type-II by construction
report("(B) rock-salt: supercells of the RELAXED RHOMBOHEDRAL AFM-II primitive (4 atoms)",
       RS_PRIM, 4, 2, [(2, 2, 2), (3, 3, 3), (4, 4, 4)])

# (A) conventional cubic, for comparison at the plan's sizes
CUB = np.eye(3) * A_RS_CUBIC
report("(A) rock-salt: CONVENTIONAL CUBIC supercells (the plan's sizes)",
       CUB, 8, 4, [(2, 2, 2), (3, 3, 3)])

report("     wurtzite: supercells of the relaxed P6_3mc cell (4 atoms)",
       WZ, 4, 2, [(2, 2, 2), (3, 3, 2), (4, 4, 3)])

print("\n" + "=" * 74)
print("READING THIS TABLE")
print("=" * 74)
print("""  Compare at EQUAL COST, not equal supercell label. The relevant question is
  how much separation each atom-count buys:

    rhombohedral (3,3,3) = 108 atoms   vs   cubic (2,2,2) = 64 atoms
    rhombohedral (2,2,2) =  32 atoms   vs   cubic (3,3,3) = 216 atoms

  A defect's interaction with its images falls off with distance, so the cell
  giving the larger minimum image distance PER ATOM is the better buy. Charge
  is not an issue here (neutral defects only) but strain and wavefunction
  overlap both scale with this distance.""")
