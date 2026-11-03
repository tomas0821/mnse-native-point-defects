#!/usr/bin/env python3
"""
Stage 3 prerequisites — three groups, one batch.

(1) BULK AFM-II REFERENCE, production settings.
    4 atoms, not 64: the AFM-II magnetic PRIMITIVE cell is rhombohedral with
    2 Mn (opposite spin) + 2 Se. Verified against the definition of type-II order
    -- 6 same-spin / 6 opposite-spin of the 12 Mn nearest neighbours. 16x cheaper
    than the conventional 2x2x2 cell, and it lets the rhombohedral distortion
    relax freely instead of being imposed on a cubic box.

(2) ORDERING HIERARCHY at the production U = 3.7574 (was checked at U_eff = 5.1).
    64-atom conventional cell, FM / AFM-I / AFM-II. Run at 70 Ry: licensed by the
    measured error cancellation (dE stable 70 -> 120 Ry to 0.02 meV/f.u.), which
    turns ~9 h into ~1.5 h.

(3) CHEMICAL-POTENTIAL REFERENCES. Needed for every formation energy; never yet
    computed. Se-rich needs only bulk Se; Mn-rich needs bulk Mn, which is what
    makes the two-limit claim (the thing distinguishing this from Aapro) possible.
    Each Mn reference is run with U AND without: applying a Hubbard U to metallic
    Mn is questionable, and the decisive test is which choice reproduces the
    experimental formation enthalpy of MnSe.

  python3 dft/gen_stage3prep.py
"""
import itertools
import math
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "stage3prep")
os.makedirs(OUT, exist_ok=True)

U_HP = 3.7574
PS = {"Mn": "Mn.pbe-spn-kjpaw_psl.0.3.1.UPF", "Se": "Se.pbe-dn-kjpaw_psl.1.0.0.UPF"}
MASS = {"Mn": 54.938044, "Se": 78.971}
A_RS = 5.5416          # our relaxed AFM-II lattice constant (70 Ry) -- starting guess


def write(tag, calc, cell, atoms, kmesh, u=U_HP, ecut=(120.0, 960.0), smear=True):
    """atoms: list of (label, x, y, z) cartesian angstrom"""
    labels = sorted({a[0] for a in atoms})
    L = ["&CONTROL", f"    calculation = '{calc}'", f"    prefix      = '{tag}'",
         "    pseudo_dir  = '/home/tomas.rojas_s/pseudo'", f"    outdir      = './tmp_{tag}'",
         "    tprnfor     = .true.", "    tstress     = .true."]
    if calc in ("relax", "vc-relax"):
        L += ["    forc_conv_thr = 1.0d-4", "    etot_conv_thr = 1.0d-5", "    nstep = 200"]
    L += ["/", "&SYSTEM", "    ibrav       = 0", f"    nat         = {len(atoms)}",
          f"    ntyp        = {len(labels)}", f"    ecutwfc     = {ecut[0]}",
          f"    ecutrho     = {ecut[1]}"]
    if smear:
        L += ["    occupations = 'smearing'", "    smearing    = 'gaussian'",
              "    degauss     = 0.01"]
    L += ["    nspin       = 2"]
    for n, lab in enumerate(labels, start=1):
        m = "+0.5" if lab.endswith("1") else ("-0.5" if lab.endswith("2") else "+0.0")
        L.append(f"    starting_magnetization({n}) = {m}")
    L += ["/", "&ELECTRONS", "    conv_thr        = 1.0d-8", "    mixing_beta     = 0.3",
          "    electron_maxstep = 200", "/"]
    if calc in ("relax", "vc-relax"):
        L += ["&IONS", "    ion_dynamics = 'bfgs'", "/"]
    if calc == "vc-relax":
        L += ["&CELL", "    cell_dynamics = 'bfgs'", "    press_conv_thr = 0.5", "/"]
    L.append("ATOMIC_SPECIES")
    for lab in labels:
        el = "Mn" if lab.startswith("Mn") else "Se"
        L.append(f"  {lab} {MASS[el]} {PS[el]}")
    if u > 0 and any(l.startswith("Mn") for l in labels):
        L.append("HUBBARD (ortho-atomic)")
        for lab in labels:
            if lab.startswith("Mn"):
                L.append(f"U {lab}-3d {u}")
    L.append("CELL_PARAMETERS angstrom")
    for v in cell:
        L.append("  {:.10f} {:.10f} {:.10f}".format(*v))
    L.append("ATOMIC_POSITIONS angstrom")
    for lab, x, y, z in atoms:
        L.append(f"  {lab} {x:.10f} {y:.10f} {z:.10f}")
    L += ["K_POINTS automatic", "  {} {} {} 0 0 0".format(*kmesh), ""]
    open(os.path.join(OUT, tag + ".in"), "w").write("\n".join(L))
    return tag


jobs = []

# ---------------- (1) AFM-II magnetic primitive cell, 4 atoms ----------------
A = np.array([[-1, 1, 0], [0, -1, 1], [1, 1, 2]]) * A_RS / 2.0
Ainv = np.linalg.inv(A)
at = []
for i, j, k in itertools.product(range(-6, 7), repeat=3):
    r = np.array([i, j, k]) * A_RS / 2.0
    f = r @ Ainv
    if all(-1e-9 <= x < 1 - 1e-9 for x in f):
        if (i + j + k) % 2 == 0:
            lab = "Mn1" if ((i + j + k) // 2) % 2 == 0 else "Mn2"
        else:
            lab = "Se"
        at.append((lab, *r))
assert len(at) == 4 and sum(a[0] == "Mn1" for a in at) == 1, at
jobs.append(write("s3_rs_AFMII_prim", "vc-relax", A, at, (8, 8, 8)))

# ---------------- (2) ordering hierarchy at the production U, 70 Ry ----------
def rocksalt64(order):
    h, out = A_RS / 2.0, []
    for i, j, k in itertools.product(range(4), repeat=3):
        if (i + j + k) % 2 == 0:
            s = (+1 if order == "FM" else (-1 if k % 2 else +1) if order == "AFMI"
                 else (-1 if ((i + j + k) // 2) % 2 else +1))
            out.append(("Mn1" if s > 0 else "Mn2", i * h, j * h, k * h))
        else:
            out.append(("Se", i * h, j * h, k * h))
    return out
cell64 = [[2 * A_RS, 0, 0], [0, 2 * A_RS, 0], [0, 0, 2 * A_RS]]
for o in ("FM", "AFMI", "AFMII"):
    jobs.append(write(f"s3_ord_{o}", "scf", cell64, rocksalt64(o), (4, 4, 4),
                      ecut=(70.0, 560.0)))

# ---------------- (3a) trigonal (grey) Se, P3_121, 3 atoms ------------------
aS, cS, xS = 4.3662, 4.9536, 0.2254
cellSe = [[aS, 0, 0], [-aS / 2, aS * math.sqrt(3) / 2, 0], [0, 0, cS]]
fr = [(xS, 0.0, 1 / 3), (0.0, xS, 2 / 3), (-xS, -xS, 0.0)]
atSe = [("Se", *(f[0] * np.array(cellSe[0]) + f[1] * np.array(cellSe[1])
                 + f[2] * np.array(cellSe[2]))) for f in fr]
jobs.append(write("s3_ref_Se", "vc-relax", cellSe, atSe, (8, 8, 6), u=0))

# ---------------- (3b) alpha-Mn, I-43m, 58 atoms ----------------------------
aM = 8.911
W = [(0.0, 0.0, 0.0), (0.3168, 0.3168, 0.3168),
     (0.3564, 0.3564, 0.0424), (0.0894, 0.0894, 0.2778)]
ops = []
for perm in itertools.permutations(range(3)):
    even = perm in [(0, 1, 2), (1, 2, 0), (2, 0, 1)]
    for s in itertools.product([1, -1], repeat=3):
        M = np.zeros((3, 3))
        for i, p in enumerate(perm):
            M[i, p] = s[i]
        d = round(np.linalg.det(M))
        if (even and d == 1) or ((not even) and d == -1):
            if not any(np.array_equal(M, N) for N in ops):
                ops.append(M)
sites = []
for w in W:
    img = set()
    for M in ops:
        for t in ([0, 0, 0], [.5, .5, .5]):
            img.add(tuple(round(x % 1.0, 6) for x in (M @ np.array(w) + np.array(t))))
    sites += sorted(img)
uniq = []
for p in sites:
    if not any(all(abs(((np.array(p) - np.array(q) + .5) % 1) - .5) < 1e-4) for q in uniq):
        uniq.append(p)
assert len(uniq) == 58, len(uniq)
# collinear AFM approximation: alternate spin by sublattice parity
atMn = [(("Mn1" if (i % 2 == 0) else "Mn2"), *(np.array(p) * aM)) for i, p in enumerate(uniq)]
cellMn = [[aM, 0, 0], [0, aM, 0], [0, 0, aM]]
jobs.append(write("s3_ref_aMn_U", "vc-relax", cellMn, atMn, (4, 4, 4), u=U_HP))
jobs.append(write("s3_ref_aMn_noU", "vc-relax", cellMn, atMn, (4, 4, 4), u=0))

# ---------------- (3c) bcc Mn AFM proxy, 2 atoms ----------------------------
aB = 2.90
cellB = [[aB, 0, 0], [0, aB, 0], [0, 0, aB]]
atB = [("Mn1", 0, 0, 0), ("Mn2", aB / 2, aB / 2, aB / 2)]
jobs.append(write("s3_ref_bccMn_U", "vc-relax", cellB, atB, (12, 12, 12), u=U_HP))
jobs.append(write("s3_ref_bccMn_noU", "vc-relax", cellB, atB, (12, 12, 12), u=0))

# cheap first, alpha-Mn last (58 atoms, metallic -- the long pole)
order = ["s3_rs_AFMII_prim", "s3_ref_Se", "s3_ref_bccMn_U", "s3_ref_bccMn_noU",
         "s3_ord_FM", "s3_ord_AFMI", "s3_ord_AFMII", "s3_ref_aMn_U", "s3_ref_aMn_noU"]

open(os.path.join(OUT, "submit_stage3prep.sh"), "w").write(f"""#!/bin/bash
#SBATCH --job-name=MnSe-s3prep
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage3prep

for tag in {' '.join(order)}; do
    if [ -f "${{tag}}.done" ]; then echo "[skip] $tag"; continue; fi
    echo "=== $tag  $(date) ==="
    mkdir -p ./tmp_${{tag}}
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${{tag}}.in > ${{tag}}.out 2>&1
    if grep -q "JOB DONE" ${{tag}}.out; then
        touch ${{tag}}.done
        echo "[ok] $tag  $(grep '!    total energy' ${{tag}}.out | tail -1)"
    else
        echo "[FAIL] $tag -- died BEFORE or AFTER '!    total energy'?"
        grep -m2 -A3 'Error in routine' ${{tag}}.out
    fi
    rm -rf ./tmp_${{tag}}
done
echo "ALL DONE $(date)"
""")
print(f"wrote {len(jobs)} inputs + submit_stage3prep.sh")
for t in order:
    n = sum(1 for _ in open(os.path.join(OUT, t + ".in")) if _.startswith("  Mn") or _.startswith("  Se"))
    print(f"   {t}")
