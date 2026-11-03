#!/bin/bash
#SBATCH --job-name=s5-fold
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=2:00:00
#SBATCH --output=slurm_fold_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/utest
# Folding emulation: reuse the ideal-u SCF density, in a COPY of the save dir so the
# existing tmp_uideal (and the bands result already in it) is left untouched.
[ -f tmp_uideal/uideal.save/charge-density.dat ] || { echo "[FAIL] no density"; exit 1; }
mkdir -p tmp_fold && cp -r tmp_uideal/uideal.save tmp_fold/
srun --mpi=pmix -n 64 pw.x -nk 8 -in uideal_fold.in > uideal_fold.out 2>&1
grep -q "JOB DONE" uideal_fold.out && echo "[ok] fold $(grep -c ' k =.*bands (ev)' uideal_fold.out) k-blocks" || echo "[FAIL] fold"
echo "FOLD DONE $(date)"
