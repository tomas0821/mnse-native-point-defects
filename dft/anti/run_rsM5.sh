#!/bin/bash
#SBATCH --job-name=rsMnSe_M5
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --output=slurm_rsMnSe_M5_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/anti
mkdir -p ./tmp_rsMnSe_M5
srun --mpi=pmix -n 64 pw.x -nk 4 -in rsMnSe_M5.in > rsMnSe_M5.out 2>&1
if grep -q "End of BFGS Geometry Optimization" rsMnSe_M5.out; then
  echo "[RELAXED] $(grep -c '^!    total energy' rsMnSe_M5.out) steps  $(grep '^!    total energy' rsMnSe_M5.out|tail -1)"
else
  echo "[INCOMPLETE] $(grep -c '^!    total energy' rsMnSe_M5.out) steps"
fi
