#!/usr/bin/env python3
"""
Generate Stage 1 Quantum ESPRESSO inputs for MnSe (rock-salt + wurtzite).

Conventions copied verbatim from the working CrN-dilution-MC/dft/gate1 setup:
  pseudo_dir, HUBBARD card style, smearing, conv_thr, mixing, SLURM header.

Structures (see 00_Master_Notebook.md Stage 1a for provenance):
  alpha-MnSe  rock-salt Fm-3m,  a = 5.4616 A      (Pradhan 2025, XRD refinement)
  w-MnSe      wurtzite  P6_3mc, a = 4.19, c = 6.68 A, u = 0.375  (Bezzerga 2025)

Magnetic orderings on the fcc Mn sublattice, in integer coords (i,j,k) of a/2
with i+j+k even:
  FM      all +
  AFM-I   sign = (-1)^k                 alternating (001) planes
  AFM-II  sign = (-1)^((i+j+k)/2)       alternating (111) planes  <- expected ground state
"""
import json
import math
import os

OUT = os.path.dirname(os.path.abspath(__file__))

# --- knobs to check against the cluster before submitting --------------------
PSEUDO_DIR = "/home/tomas.rojas_s/pseudo"
PS_MN = "Mn.pbe-spn-kjpaw_psl.0.3.1.UPF"
PS_SE = "Se.pbe-dn-kjpaw_psl.1.0.0.UPF"
# Provisional U for this batch. NOT arbitrary and NOT a claim:
# the HUBBARD card below specifies U only, which in QE is the Dudarev-style
# effective U_eff = U - J. The anchor (Grzybowski 2024) uses the Liechtenstein
# scheme with U = 6.0 eV and J_H = 0.15U = 0.9 eV, whose Dudarev equivalent is
# U_eff = 6.0 - 0.9 = 5.1 eV. So 5.1 is the value that most nearly mimics the
# anchor under the syntax proven to run on this cluster -- which makes Batch B
# directly comparable to the anchor's AFM-vs-FM result (its Fig. 6d).
# (It coincides with Bezzerga's U = 5.1 eV; that appears to be a coincidence.)
U_PROV = 5.1
A_RS = 5.4616           # alpha-MnSe cubic lattice constant, A
A_WZ, C_WZ, U_WZ = 4.19, 6.68, 0.375

MASS = {"Mn": 54.938044, "Se": 78.971}


# --- structure builders -------------------------------------------------------
def rocksalt(n=1):
    """n x n x n conventional cubic cells. Returns (cell, [(elem,i,j,k)]) with
    i,j,k integer coords in units of a/2."""
    L = 2 * n
    mn, se = [], []
    for i in range(L):
        for j in range(L):
            for k in range(L):
                (mn if (i + j + k) % 2 == 0 else se).append((i, j, k))
    cell = [[A_RS * n, 0, 0], [0, A_RS * n, 0], [0, 0, A_RS * n]]
    return cell, mn, se


def spin_sign(order, i, j, k):
    if order == "FM":
        return +1
    if order == "AFMI":
        return -1 if k % 2 else +1
    if order == "AFMII":
        return -1 if ((i + j + k) // 2) % 2 else +1
    raise ValueError(order)


def wurtzite():
    cell = [[A_WZ, 0.0, 0.0],
            [-A_WZ / 2, A_WZ * math.sqrt(3) / 2, 0.0],
            [0.0, 0.0, C_WZ]]
    mn = [(1 / 3, 2 / 3, 0.0), (2 / 3, 1 / 3, 0.5)]
    se = [(1 / 3, 2 / 3, U_WZ), (2 / 3, 1 / 3, U_WZ + 0.5)]
    return cell, mn, se


# --- input writer -------------------------------------------------------------
def write_pw(tag, calc, cell, species_pos, kmesh, ecutwfc, ecutrho, u=U_PROV,
             two_mn=True, press_conv=0.5):
    """species_pos: list of (label, x, y, z) in ANGSTROM (cartesian)."""
    labels = sorted({s[0] for s in species_pos})
    ntyp = len(labels)
    lines = []
    lines.append("&CONTROL")
    lines.append(f"    calculation = '{calc}'")
    lines.append(f"    prefix      = '{tag}'")
    lines.append(f"    pseudo_dir  = '{PSEUDO_DIR}'")
    lines.append("    outdir      = './tmp'")
    lines.append("    tprnfor     = .true.")
    lines.append("    tstress     = .true.")
    if calc in ("relax", "vc-relax"):
        lines.append("    forc_conv_thr = 1.0d-4")
        lines.append("    etot_conv_thr = 1.0d-5")
        lines.append("    nstep       = 200")
    lines.append("/")
    lines.append("&SYSTEM")
    lines.append("    ibrav       = 0")
    lines.append(f"    nat         = {len(species_pos)}")
    lines.append(f"    ntyp        = {ntyp}")
    lines.append(f"    ecutwfc     = {ecutwfc:.1f}")
    lines.append(f"    ecutrho     = {ecutrho:.1f}")
    lines.append("    occupations = 'smearing'")
    lines.append("    smearing    = 'gaussian'")
    lines.append("    degauss     = 0.01")
    lines.append("    nspin       = 2")
    for n, lab in enumerate(labels, start=1):
        if lab == "Mn1":
            lines.append(f"    starting_magnetization({n}) = +0.5")
        elif lab == "Mn2":
            lines.append(f"    starting_magnetization({n}) = -0.5")
        else:
            lines.append(f"    starting_magnetization({n}) = +0.0")
    lines.append("/")
    lines.append("&ELECTRONS")
    lines.append("    conv_thr        = 1.0d-8")
    lines.append("    mixing_beta     = 0.3")
    lines.append("    electron_maxstep = 200")
    lines.append("/")
    if calc in ("relax", "vc-relax"):
        lines.append("&IONS")
        lines.append("    ion_dynamics = 'bfgs'")
        lines.append("/")
    if calc == "vc-relax":
        lines.append("&CELL")
        lines.append("    cell_dynamics = 'bfgs'")
        lines.append(f"    press_conv_thr = {press_conv}")
        lines.append("/")
    lines.append("ATOMIC_SPECIES")
    for lab in labels:
        el = "Mn" if lab.startswith("Mn") else "Se"
        lines.append(f"  {lab} {MASS[el]} {PS_MN if el=='Mn' else PS_SE}")
    if u > 0:
        lines.append("HUBBARD (ortho-atomic)")
        for lab in labels:
            if lab.startswith("Mn"):
                lines.append(f"U {lab}-3d {u}")
    lines.append("CELL_PARAMETERS angstrom")
    for v in cell:
        lines.append("  {:.10f} {:.10f} {:.10f}".format(*v))
    lines.append("ATOMIC_POSITIONS angstrom")
    for lab, x, y, z in species_pos:
        lines.append(f"  {lab} {x:.10f} {y:.10f} {z:.10f}")
    lines.append("K_POINTS automatic")
    lines.append("  {} {} {} 0 0 0".format(*kmesh))
    lines.append("")
    path = os.path.join(OUT, "stage1", f"{tag}.in")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(lines))
    return f"{tag}.in"


def rs_positions(n, order):
    cell, mn, se = rocksalt(n)
    h = A_RS / 2
    pos = []
    for (i, j, k) in mn:
        lab = "Mn1" if spin_sign(order, i, j, k) > 0 else "Mn2"
        pos.append((lab, i * h, j * h, k * h))
    for (i, j, k) in se:
        pos.append(("Se", i * h, j * h, k * h))
    return cell, pos


def wz_positions(order):
    cell, mn, se = wurtzite()
    def cart(f):
        return [sum(f[i] * cell[i][d] for i in range(3)) for d in range(3)]
    pos = []
    for n, f in enumerate(mn):
        lab = "Mn1" if (order == "FM" or n == 0) else "Mn2"
        x, y, z = cart(f)
        pos.append((lab, x, y, z))
    for f in se:
        x, y, z = cart(f)
        pos.append(("Se", x, y, z))
    return cell, pos


# --- build the batches --------------------------------------------------------
jobs = {"A_convergence": [], "B_ordering": [], "C_relax": []}

# Batch A1: cutoff scan, small cells, fixed dense-ish k
# NOTE: uses AFM-I, not AFM-II. Type-II order does NOT fit the 8-atom conventional
# cubic cell -- its (111) alternation needs the cell doubled to 2a (the 64-atom cell),
# and forcing it here gives 1 up / 3 down, i.e. a net moment. AFM-I fits exactly
# (2 up / 2 down) and cutoff convergence is insensitive to which AFM pattern is used.
cell, pos = rs_positions(1, "AFMI")           # 8-atom conventional cubic
for ec in (40, 50, 60, 70, 80, 90):
    jobs["A_convergence"].append(
        write_pw(f"convA_rs_ecut{ec}", "scf", cell, pos, (6, 6, 6), ec, 8 * ec))
cellw, posw = wz_positions("AFM")             # 4-atom wurtzite primitive
for ec in (40, 50, 60, 70, 80, 90):
    jobs["A_convergence"].append(
        write_pw(f"convA_wz_ecut{ec}", "scf", cellw, posw, (8, 8, 6), ec, 8 * ec))

# Batch A2: k-mesh scan at a provisional cutoff (70 Ry) -- revise after A1
for k in (4, 6, 8, 10):
    jobs["A_convergence"].append(
        write_pw(f"convA_rs_k{k}", "scf", cell, pos, (k, k, k), 70, 560))
for k, kz in ((6, 4), (8, 6), (10, 6), (12, 8)):
    jobs["A_convergence"].append(
        write_pw(f"convA_wz_k{k}", "scf", cellw, posw, (k, k, kz), 70, 560))

# Batch B: magnetic ordering comparison at the experimental geometry
for order in ("FM", "AFMI", "AFMII"):
    c, p = rs_positions(2, order)             # 64-atom 2x2x2 conventional
    jobs["B_ordering"].append(
        write_pw(f"ordB_rs_{order}", "scf", c, p, (4, 4, 4), 70, 560))
for order in ("FM", "AFM"):
    c, p = wz_positions(order)
    jobs["B_ordering"].append(
        write_pw(f"ordB_wz_{order}", "scf", c, p, (10, 10, 6), 70, 560))

# Batch C: variable-cell relaxation in the expected ground states.
# Cheap wurtzite job FIRST so it is guaranteed to finish overnight; the 64-atom
# rock-salt AFM-II vc-relax goes LAST because it is by far the most expensive and
# everything else should complete even if it does not.
# AFM-II must use the 64-atom cell (see note in Batch A1). ibrav=0 + vc-relax lets
# the cell distort rhombohedrally -- the explicit test of magnetostructural coupling.
c, p = wz_positions("AFM")
jobs["C_relax"].append(write_pw("relC_wz_AFM", "vc-relax", c, p, (10, 10, 6), 70, 560))
c, p = rs_positions(1, "AFMI")
jobs["C_relax"].append(write_pw("relC_rs_AFMI_quick", "vc-relax", c, p, (6, 6, 6), 70, 560))
c, p = rs_positions(2, "AFMII")
jobs["C_relax"].append(write_pw("relC_rs_AFMII", "vc-relax", c, p, (4, 4, 4), 70, 560))

order_all = jobs["A_convergence"] + jobs["B_ordering"] + jobs["C_relax"]
tags = [f[:-3] for f in order_all]

with open(os.path.join(OUT, "stage1", "manifest.json"), "w") as f:
    json.dump({"batches": {k: [t[:-3] for t in v] for k, v in jobs.items()},
               "U_provisional": U_PROV,
               "a_rocksalt": A_RS, "a_wurtzite": A_WZ, "c_wurtzite": C_WZ,
               "pseudo": {"Mn": PS_MN, "Se": PS_SE}}, f, indent=2)

SUBMIT = f"""#!/bin/bash
#SBATCH --job-name=MnSe-stage1
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

cd /scratch/global/tomas.rojas_s/mnse/stage1

for tag in {' '.join(tags)}; do
    if [ -f "${{tag}}.done" ]; then
        echo "[skip] $tag already done"
        continue
    fi
    echo "=== $tag  $(date) ==="
    # mkdir BEFORE the run. The loop deletes ./tmp at the end of each iteration,
    # and 64 MPI ranks then race to recreate it -- some lose and pw.x dies with
    # "unable to create directory ./tmp/<prefix>.save/" AFTER the SCF has already
    # converged. Cost us a 4 h 64-atom job on 2026-08-07. Creating it up front,
    # serially, removes the race.
    mkdir -p ./tmp
    # This cluster PROHIBITS mpirun/mpiexec (SLURMMASTER rejects the job at
    # submit time). srun --mpi=pmix is mandatory. The CrN gate1 script used
    # mpirun and predates this policy -- do not copy it verbatim again.
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${{tag}}.in > ${{tag}}.out 2>&1
    if grep -q "JOB DONE" ${{tag}}.out; then
        touch ${{tag}}.done
        echo "[ok] $tag  $(grep '!    total energy' ${{tag}}.out | tail -1)"
    else
        echo "[FAIL] $tag -- see ${{tag}}.out"
    fi
    rm -rf ./tmp
done
echo "ALL DONE $(date)"
"""
with open(os.path.join(OUT, "stage1", "submit_stage1.sh"), "w") as f:
    f.write(SUBMIT)

print(f"wrote {len(order_all)} inputs + submit_stage1.sh + manifest.json")
for b, v in jobs.items():
    print(f"  {b}: {len(v)}")
