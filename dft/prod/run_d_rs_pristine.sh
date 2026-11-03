#!/bin/bash
#SBATCH --job-name=d_rs_pristine
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=08:00:00
#SBATCH --output=slurm_d_rs_pristine_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/prod
mkdir -p ./tmp_d_rs_pristine_s120
srun --mpi=pmix -n 64 pw.x -nk 4 -in d_rs_pristine_s120.in > d_rs_pristine_s120.out 2>&1
grep -q "JOB DONE" d_rs_pristine_s120.out && ! grep -q "convergence NOT achieved" d_rs_pristine_s120.out && { touch d_rs_pristine.done; echo "[ok] d_rs_pristine  $(grep '^!    total energy' d_rs_pristine_s120.out|tail -1)"; } || echo "[FAIL] d_rs_pristine"
rm -rf ./tmp_d_rs_pristine_s120
