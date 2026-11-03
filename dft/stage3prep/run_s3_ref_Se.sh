#!/bin/bash
#SBATCH --job-name=s3_ref_Se
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm_s3_ref_Se_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage3prep

if [ -f "s3_ref_Se.done" ]; then echo "[skip] s3_ref_Se already complete"; exit 0; fi
mkdir -p ./tmp_s3_ref_Se
srun --mpi=pmix -n 8 pw.x -nk 8 -in s3_ref_Se.in > s3_ref_Se.out 2>&1
if grep -q "JOB DONE" s3_ref_Se.out; then
    touch s3_ref_Se.done
    echo "[ok] s3_ref_Se  $(grep '!    total energy' s3_ref_Se.out | tail -1)"
else
    echo "[FAIL] s3_ref_Se -- died BEFORE or AFTER '!    total energy'?"
    grep -m2 -A3 'Error in routine' s3_ref_Se.out
fi
rm -rf ./tmp_s3_ref_Se
