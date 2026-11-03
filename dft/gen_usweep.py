#!/usr/bin/env python3
"""
Stage 2 — the U sweep.

PURPOSE: diagnose why our relaxed volumes sit 2-4% above the VASP literature
(a = 5.5791 A vs 5.4616 experimental; w-MnSe a = 4.2937 vs 4.19 Bezzerga / 4.12
anchor), while c/a agrees. Suspicion: the Hubbard U value and/or scheme.

This is NOT hp.x. hp.x *determines* U ab initio; it cannot tell you why you
disagree with someone else's VASP number, and it must be evaluated at a settled
geometry. Sweep first, hp.x after. See the notebook.

Built from the converged Stage 1 settings: ecutwfc 120 / ecutrho 960, k as
established. All jobs are vc-relax, because the lattice constant is the observable
under investigation.

  A_dudarev   : U_eff = 0, 2, 3, 4, 5.1, 6      (both polymorphs)   12 jobs
  B_scheme    : Liechtenstein (U=6.0, J=0.9)    (both polymorphs)    2 jobs
                -- the anchor's EXACT parameters. Paired against the existing
                Dudarev U_eff = 5.1 run: since 6.0 - 0.9 = 5.1, any difference
                between them is the SCHEME, not the value. That is the controlled
                experiment.

  python3 dft/gen_usweep.py
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "stage1")
OUT = os.path.join(HERE, "usweep")
os.makedirs(OUT, exist_ok=True)

TEMPLATES = {"rs": "relC_rs_AFMI_quick_c120.in", "wz": "relC_wz_AFM_c120.in"}
U_DUDAREV = [0.0, 2.0, 3.0, 4.0, 5.1, 6.0]
LIECHT = (6.0, 0.9)


def rewrite(src_text, tag, u, j=None):
    t = src_text
    t = re.sub(r"prefix\s*=\s*'[^']*'", f"prefix      = '{tag}'", t)
    # strip the existing HUBBARD block (card + its U lines), rebuild it
    lines, out, in_hub = t.splitlines(), [], False
    for ln in lines:
        if ln.startswith("HUBBARD"):
            in_hub = True
            continue
        if in_hub:
            if re.match(r"^[UJ]\s+Mn\d?-3d", ln):
                continue
            in_hub = False
        out.append(ln)
    t = "\n".join(out)
    if u > 0:
        block = ["HUBBARD (ortho-atomic)"]
        for sp in ("Mn1", "Mn2"):
            if f"  {sp} " in t:
                block.append(f"U {sp}-3d {u}")
                if j is not None:
                    block.append(f"J {sp}-3d {j}")
        t = t.replace("CELL_PARAMETERS", "\n".join(block) + "\nCELL_PARAMETERS")
    return t.rstrip() + "\n"


jobs = []
for sysname, fn in TEMPLATES.items():
    src = open(os.path.join(SRC, fn)).read()
    for u in U_DUDAREV:
        tag = f"uA_{sysname}_U{str(u).replace('.', 'p')}"
        open(os.path.join(OUT, tag + ".in"), "w").write(rewrite(src, tag, u))
        jobs.append(tag)
    u, j = LIECHT
    tag = f"uB_{sysname}_L_U{str(u).replace('.', 'p')}_J{str(j).replace('.', 'p')}"
    open(os.path.join(OUT, tag + ".in"), "w").write(rewrite(src, tag, u, j))
    jobs.append(tag)

# cheap jobs first: wurtzite (4 atoms) before rock-salt (8 atoms)
jobs.sort(key=lambda t: (0 if "_wz_" in t else 1, t))

SUBMIT = f"""#!/bin/bash
#SBATCH --job-name=MnSe-usweep
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu

cd /scratch/global/tomas.rojas_s/mnse/usweep

for tag in {' '.join(jobs)}; do
    if [ -f "${{tag}}.done" ]; then echo "[skip] $tag"; continue; fi
    echo "=== $tag  $(date) ==="
    mkdir -p ./tmp          # serial mkdir: avoids the 64-rank race that killed
                            # ordB_rs_FM_c120 after its SCF had already converged
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${{tag}}.in > ${{tag}}.out 2>&1
    if grep -q "JOB DONE" ${{tag}}.out; then
        touch ${{tag}}.done
        echo "[ok] $tag  $(grep '!    total energy' ${{tag}}.out | tail -1)"
    else
        echo "[FAIL] $tag -- check whether it died BEFORE or AFTER '!    total energy'"
    fi
    rm -rf ./tmp
done
echo "ALL DONE $(date)"
"""
open(os.path.join(OUT, "submit_usweep.sh"), "w").write(SUBMIT)
print(f"wrote {len(jobs)} inputs + submit_usweep.sh")
for t in jobs:
    print("  ", t)
