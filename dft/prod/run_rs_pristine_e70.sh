#!/bin/bash
#SBATCH --job-name=rs64-e70
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=6:00:00
#SBATCH --output=slurm_rs64_e70_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/prod
# Finite-size reference: 64-atom pristine alpha-MnSe at 70 Ry in the PRODUCTION cell, matched to
# the 70 Ry defect relaxations (d_rs_*_r70) and to the 128-atom fs128 runs. Scratch deleted after.
mkdir -p ./tmp_d_rs_pristine_e70
srun --mpi=pmix -n 32 pw.x -nk 4 -in d_rs_pristine_e70.in > d_rs_pristine_e70.out 2>&1
grep -q "JOB DONE" d_rs_pristine_e70.out && echo "[ok] E=$(grep '^!' d_rs_pristine_e70.out | tail -1 | awk '{print $5}')" || echo "[FAIL]"
rm -rf ./tmp_d_rs_pristine_e70
