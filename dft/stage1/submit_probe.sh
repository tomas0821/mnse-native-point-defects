#!/bin/bash
#SBATCH --job-name=MnSe-probe
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:20:00
#SBATCH --output=slurm_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage1
mkdir -p ./tmp_probe
srun --mpi=pmix -n 16 pw.x -in probe_UJ.in > probe_UJ.out 2>&1
echo "exit=$?"
grep -iE "Hubbard|lda_plus_u|Error|%%%%" probe_UJ.out | head -12
