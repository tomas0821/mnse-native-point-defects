#!/bin/bash
#SBATCH --job-name=MnSe-wzb
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=00:30:00
#SBATCH --output=slurm_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/bands
cat b_wz_bands.head geom_b_wz.txt b_wz_kpath.txt > b_wz_bands.in
srun --mpi=pmix -n 16 pw.x -nk 4 -in b_wz_bands.in > b_wz_bands.out 2>&1
grep -q "JOB DONE" b_wz_bands.out && echo "bands OK, eigenvalue blocks: $(grep -c 'bands (ev)' b_wz_bands.out)" || echo FAILED
