#!/bin/bash
#SBATCH --job-name=aMn2_NM_noU
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=12:00:00
#SBATCH --output=slurm_aMn2_NM_noU_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage3prep
mkdir -p ./tmp_aMn2_NM_noU
srun --mpi=pmix -n 64 pw.x -nk 4 -in aMn2_NM_noU.in > aMn2_NM_noU.out 2>&1
grep -q "JOB DONE" aMn2_NM_noU.out && { touch aMn2_NM_noU.done; echo "[ok] aMn2_NM_noU $(grep '!    total energy' aMn2_NM_noU.out|tail -1)"; }                           || echo "[FAIL] aMn2_NM_noU $(grep -c 'iteration #' aMn2_NM_noU.out) iterations, no convergence"
rm -rf ./tmp_aMn2_NM_noU
