#!/bin/bash
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
