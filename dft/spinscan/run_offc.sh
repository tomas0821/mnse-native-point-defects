#!/bin/bash
#SBATCH --job-name=rsSeMn-offc
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=36:00:00
#SBATCH --output=slurm_%x_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
export OMP_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
srun --mpi=pmix -n 64 pw.x -nk 4 -in rsSeMn_offc.in > rsSeMn_offc.out 2>&1
rm -rf tmp_rsSeMn_offc
