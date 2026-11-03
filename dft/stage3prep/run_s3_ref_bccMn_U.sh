#!/bin/bash
#SBATCH --job-name=s3_ref_bccMn_U
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm_s3_ref_bccMn_U_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage3prep

if [ -f "s3_ref_bccMn_U.done" ]; then echo "[skip] s3_ref_bccMn_U already complete"; exit 0; fi
mkdir -p ./tmp_s3_ref_bccMn_U
srun --mpi=pmix -n 8 pw.x -nk 8 -in s3_ref_bccMn_U.in > s3_ref_bccMn_U.out 2>&1
if grep -q "JOB DONE" s3_ref_bccMn_U.out; then
    touch s3_ref_bccMn_U.done
    echo "[ok] s3_ref_bccMn_U  $(grep '!    total energy' s3_ref_bccMn_U.out | tail -1)"
else
    echo "[FAIL] s3_ref_bccMn_U -- died BEFORE or AFTER '!    total energy'?"
    grep -m2 -A3 'Error in routine' s3_ref_bccMn_U.out
fi
rm -rf ./tmp_s3_ref_bccMn_U
