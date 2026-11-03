#!/usr/bin/env python3
"""
Stage 1b — the two gates that have never actually been measured:
  (i)  the band gaps, and specifically the GAP DIFFERENCE (target ~1.5 eV, the
       anchor is explicit that the difference is trustworthy and the absolutes
       are not), and
  (ii) the altermagnetic spin splitting of w-MnSe along L1-Gamma-L2
       (target ~245 meV, Bezzerga; the anchor publishes no number).

Everything runs at the ab initio U = 3.76 eV from hp.x.

Three chained steps per polymorph, in one job:
  relax  vc-relax at U=3.76, smearing        -> production geometry
  scf    fixed occupations at that geometry  -> QE then PRINTS the HOMO and LUMO
                                                directly, which is the cleanest
                                                way to get a gap (no band path,
                                                no folding ambiguity)
  bands  along L1-Gamma-L2 (wurtzite)        -> the spin splitting

The geometry is spliced from the vc-relax output at runtime ("Begin final
coordinates" block), so the later steps always use the relaxed cell.

  python3 dft/gen_bands.py
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "stage1")
OUT = os.path.join(HERE, "bands")
os.makedirs(OUT, exist_ok=True)

U_HP = 3.7574

# wurtzite: L1 -> Gamma -> L2. The six L points of the hexagonal BZ split into
# two sets related by the altermagnetic operation (6-fold rotation about c plus a
# non-primitive translation), so L1 and L2 must show splitting of OPPOSITE sign.
# Taking them as the two path endpoints means one run gives both.
KPATH_WZ = """K_POINTS crystal_b
3
  0.5000 0.0000 0.5000 60   ! L1
  0.0000 0.0000 0.0000 60   ! Gamma
  0.0000 0.5000 0.5000 1    ! L2
"""
# rock-salt (8-atom conventional cubic cell -> simple-cubic BZ). Used only to
# locate the gap; bands here are folded, so do NOT read dispersion from it.
KPATH_RS = """K_POINTS crystal_b
4
  0.0000 0.0000 0.0000 40   ! Gamma
  0.5000 0.0000 0.0000 40   ! X
  0.5000 0.5000 0.0000 40   ! M
  0.0000 0.0000 0.0000 1    ! Gamma
"""

SYSTEMS = [("relC_wz_AFM_c120.in", "b_wz", 40, KPATH_WZ),
           ("relC_rs_AFMI_quick_c120.in", "b_rs", 75, KPATH_RS)]


def strip_blocks(t):
    t = re.sub(r"&IONS.*?\n/\n", "", t, flags=re.S)
    t = re.sub(r"&CELL.*?\n/\n", "", t, flags=re.S)
    for k in ("forc_conv_thr", "etot_conv_thr", "nstep"):
        t = re.sub(rf"\s*{k}\s*=.*\n", "\n", t)
    return t


def set_u(t, u):
    t = re.sub(r"^U (Mn\d?)-3d .*$", lambda m: f"U {m.group(1)}-3d {u}", t, flags=re.M)
    return t


def header_only(t):
    """everything before CELL_PARAMETERS -- geometry gets appended at runtime"""
    return t.split("CELL_PARAMETERS")[0]


for fn, tag, nbnd, kpath in SYSTEMS:
    src = set_u(open(os.path.join(SRC, fn)).read(), U_HP)

    # --- step 1: vc-relax at the ab initio U (smearing: robust for relaxation) ---
    t = re.sub(r"prefix\s*=\s*'[^']*'", f"prefix      = '{tag}'", src)
    t = re.sub(r"outdir\s*=\s*'[^']*'", f"outdir      = './tmp_{tag}'", t)
    open(os.path.join(OUT, f"{tag}_relax.in"), "w").write(t)

    # --- step 2: fixed-occupation scf -> QE prints highest occupied / lowest unoccupied
    t2 = strip_blocks(t)
    t2 = re.sub(r"calculation\s*=\s*'[^']*'", "calculation = 'scf'", t2)
    t2 = t2.replace("    occupations = 'smearing'\n", "    occupations = 'fixed'\n")
    t2 = re.sub(r"\s*smearing\s*=.*\n", "\n", t2)
    t2 = re.sub(r"\s*degauss\s*=.*\n", "\n", t2)
    t2 = t2.replace("    nspin       = 2",
                    f"    nbnd        = {nbnd}\n    tot_magnetization = 0.0\n    nspin       = 2")
    open(os.path.join(OUT, f"{tag}_scf.head"), "w").write(header_only(t2))

    # --- step 3: bands along the path ---
    t3 = re.sub(r"calculation\s*=\s*'[^']*'", "calculation = 'bands'", t2)
    open(os.path.join(OUT, f"{tag}_bands.head"), "w").write(header_only(t3))
    open(os.path.join(OUT, f"{tag}_kpath.txt"), "w").write(kpath)

SUBMIT = """#!/bin/bash
#SBATCH --job-name=MnSe-bands
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=slurm_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/bands

# splice the relaxed cell+positions out of a vc-relax output
extract_geom () {   # $1 = vc-relax .out
    awk '/Begin final coordinates/,/End final coordinates/' "$1" \\
      | sed -n '/CELL_PARAMETERS/,$p' \\
      | grep -v "End final coordinates"
}

for tag in b_wz b_rs; do
    echo "######## $tag  $(date) ########"
    mkdir -p ./tmp_${tag}

    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}_relax.in > ${tag}_relax.out 2>&1
    if ! grep -q "JOB DONE" ${tag}_relax.out; then echo "  [FAIL] relax"; continue; fi
    echo "  relax: $(grep '^!    total energy' ${tag}_relax.out | tail -1)"

    extract_geom ${tag}_relax.out > geom_${tag}.txt
    if [ ! -s geom_${tag}.txt ]; then echo "  [FAIL] no final geometry"; continue; fi

    cat ${tag}_scf.head   geom_${tag}.txt > ${tag}_scf.in
    printf 'K_POINTS automatic\\n' >> ${tag}_scf.in
    grep -A1 "K_POINTS automatic" ${tag}_relax.in | tail -1 >> ${tag}_scf.in

    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}_scf.in > ${tag}_scf.out 2>&1
    if grep -q "JOB DONE" ${tag}_scf.out; then
        echo "  GAP:"; grep -iE "highest occupied|lowest unoccupied" ${tag}_scf.out | tail -2
    else
        echo "  [FAIL] scf"; grep -A3 'Error in routine' ${tag}_scf.out | head -6; continue
    fi

    cat ${tag}_bands.head geom_${tag}.txt ${tag}_kpath.txt > ${tag}_bands.in
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}_bands.in > ${tag}_bands.out 2>&1
    grep -q "JOB DONE" ${tag}_bands.out && echo "  [ok] bands done" || echo "  [FAIL] bands"
done
echo "ALL DONE $(date)"
"""
open(os.path.join(OUT, "submit_bands.sh"), "w").write(SUBMIT)
print("wrote inputs + submit_bands.sh in", OUT)
for f in sorted(os.listdir(OUT)):
    print("  ", f)
