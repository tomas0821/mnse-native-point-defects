#!/bin/bash
# Unrelaxed 72-atom wurtzite defect cell (ideal positions of the production relaxation input), 70/560 Ry.
# Usage: sbatch --job-name=ur-TAG run_ur.sh NAME
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%x_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
export OMP_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
srun --mpi=pmix -n 64 pw.x -nk 4 -in "$1.in" > "$1.out" 2>&1
echo "$1: $(grep -c 'iteration #' "$1.out") iterations, $(grep '^!' "$1.out" | tail -1), $(grep 'total magnetization' "$1.out" | tail -1)"
rm -rf "tmp_$1"
