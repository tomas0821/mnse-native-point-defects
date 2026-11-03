#!/usr/bin/env python3
"""
Fixed-occupation template + hp.x linear-response U.

WHY: every input so far uses Gaussian smearing, inherited from the CrN template
where it is right (CrN is nearly metallic). MnSe is a magnetic SEMICONDUCTOR --
hp.x refused to run and told us so:
    DOS(E_Fermi) = 0.1256E-37 ... the system has a gap and should NOT be treated
    as a metal. Solution (for magnetic insulators): use the 2-step scf procedure.

So each system gets three inputs:
  step1  scf, smearing            -> converges the density the easy way
  step2  scf, occupations='fixed' -> restarts from step1's density (startingpot='file')
                                     with tot_magnetization pinned. This is the
                                     state hp.x needs, and its energy must match
                                     step1 to within noise or something is wrong.
  hp     hp.x linear response     -> the actual U

U is a local atomic property, so it is determined on the CHEAP cells (4-atom
wurtzite, 8-atom rock-salt AFM-I) rather than the 64-atom AFM-II ground-state
cell. If the two agree, U is transferable and that is itself worth reporting.

Valence from the pseudopotentials actually in use: Mn 15 e-, Se 16 e-.

  python3 dft/gen_hp.py
"""
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "stage1")
OUT = os.path.join(HERE, "hp")
os.makedirs(OUT, exist_ok=True)

ZVAL = {"Mn": 15, "Se": 16}
# (template, tag, nbnd, q-mesh)  -- q must divide the k-mesh in the template
SYSTEMS = [
    ("ordB_wz_AFM_c120.in", "hp_wz", 40, (2, 2, 2)),
    ("relC_rs_AFMI_quick_c120.in", "hp_rs", 75, (2, 2, 2)),
]


def count_electrons(text):
    body = text.split("ATOMIC_POSITIONS")[1].split("K_POINTS")[0]
    n = 0
    for ln in body.splitlines():
        p = ln.split()
        if len(p) == 4:
            n += ZVAL["Mn" if p[0].startswith("Mn") else "Se"]
    return n


def make(src_text, tag, nbnd, step):
    t = src_text
    # every step of a system shares prefix+outdir so step2 and hp.x can restart
    t = re.sub(r"prefix\s*=\s*'[^']*'", f"prefix      = '{tag}'", t)
    t = re.sub(r"outdir\s*=\s*'[^']*'", f"outdir      = './tmp_{tag}'", t)
    t = re.sub(r"calculation\s*=\s*'[^']*'", "calculation = 'scf'", t)
    # a 'relax'/'vc-relax' template carries &IONS/&CELL; for scf they are stray
    # cards -- leaving them in is what broke the earlier +U+V probe
    t = re.sub(r"&IONS.*?\n/\n", "", t, flags=re.S)
    t = re.sub(r"&CELL.*?\n/\n", "", t, flags=re.S)
    for k in ("forc_conv_thr", "etot_conv_thr", "nstep"):
        t = re.sub(rf"\s*{k}\s*=.*\n", "\n", t)
    # hp.x needs the wavefunctions on disk
    t = t.replace("&CONTROL", "&CONTROL\n    disk_io     = 'medium'")
    # nbnd must be set in BOTH steps and be identical. Step 1 left to itself picks
    # the smeared-metal default (37 for the 4-atom cell, 74 for the 8-atom), and
    # step 2 then aborts with "The number of bands for this run is 40, but only 37
    # were read from file".
    t = t.replace("    nspin       = 2", f"    nbnd        = {nbnd}\n    nspin       = 2")
    if step == 2:
        t = t.replace("    occupations = 'smearing'\n", "    occupations = 'fixed'\n")
        t = re.sub(r"\s*smearing\s*=.*\n", "\n", t)
        t = re.sub(r"\s*degauss\s*=.*\n", "\n", t)
        nelec = count_electrons(t)
        t = t.replace("    nspin       = 2",
                      f"    tot_magnetization = 0.0   ! {nelec} e-, compensated AFM\n"
                      "    nspin       = 2")
        # restart from the DENSITY only. Reusing the wavefunctions as well adds a
        # second way to trip the same band-count mismatch for no benefit.
        t = t.replace("&ELECTRONS", "&ELECTRONS\n    startingpot = 'file'")
    return t


made = []
for fn, tag, nbnd, nq in SYSTEMS:
    src = open(os.path.join(SRC, fn)).read()
    for step in (1, 2):
        name = f"{tag}_step{step}"
        open(os.path.join(OUT, name + ".in"), "w").write(make(src, tag, nbnd, step))
        made.append(name)
    open(os.path.join(OUT, tag + "_hp.in"), "w").write(
        f"""&INPUTHP
  prefix       = '{tag}'
  outdir       = './tmp_{tag}'
  nq1 = {nq[0]}, nq2 = {nq[1]}, nq3 = {nq[2]}
  conv_thr_chi = 1.0d-5
  iverbosity   = 1
/
""")
    made.append(tag + "_hp")

SUBMIT = """#!/bin/bash
#SBATCH --job-name=MnSe-hp
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
cd /scratch/global/tomas.rojas_s/mnse/hp

for tag in hp_wz hp_rs; do
    echo "######## $tag  $(date) ########"
    mkdir -p ./tmp_${tag}

    # step 1 -- smeared scf, gets the density cheaply
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}_step1.in > ${tag}_step1.out 2>&1
    E1=$(grep '^!    total energy' ${tag}_step1.out | tail -1)
    echo "  step1 (smearing): $E1"

    # step 2 -- fixed occupations, restarts from step1. This is what hp.x needs.
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}_step2.in > ${tag}_step2.out 2>&1
    E2=$(grep '^!    total energy' ${tag}_step2.out | tail -1)
    echo "  step2 (fixed)   : $E2"
    if ! grep -q "JOB DONE" ${tag}_step2.out; then
        echo "  [FAIL] step2 -- hp.x cannot run, skipping"
        grep -A3 'Error in routine' ${tag}_step2.out | head -6
        continue
    fi

    # hp.x -- the actual linear-response U
    srun --mpi=pmix -n 64 hp.x -in ${tag}_hp.in > ${tag}_hp.out 2>&1
    if grep -q "JOB DONE" ${tag}_hp.out && ! grep -qi "Stopping" ${tag}_hp.out; then
        echo "  [ok] hp.x done"
        grep -A12 -iE "Hubbard U parameters|Computed Hubbard parameters" ${tag}_hp.out | head -16
    else
        echo "  [FAIL] hp.x"
        grep -iE "WARNING|Stopping|Error" ${tag}_hp.out | head -8
    fi
    # keep tmp_* -- hp.x results live there; do NOT rm
done
echo "ALL DONE $(date)"
"""
open(os.path.join(OUT, "submit_hp.sh"), "w").write(SUBMIT)
print(f"wrote {len(made)} inputs + submit_hp.sh")
for m in made:
    print("  ", m)
