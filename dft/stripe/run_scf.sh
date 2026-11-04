#!/bin/bash
#SBATCH --job-name=stripe-scf
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --output=slurm_%x_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
export OMP_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
for n in AFMsa Atype AFMsb FM; do
  srun --mpi=pmix -n 32 pw.x -nk 4 -in $n.in > $n.out 2>&1
  rm -rf tmp_$n
done
