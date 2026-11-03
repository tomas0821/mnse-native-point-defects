#!/bin/bash
#SBATCH --job-name=s3_ref_aMn_noU
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=20:00:00
#SBATCH --output=slurm_s3_ref_aMn_noU_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage3prep

if [ -f "s3_ref_aMn_noU.done" ]; then echo "[skip] s3_ref_aMn_noU already complete"; exit 0; fi
mkdir -p ./tmp_s3_ref_aMn_noU
srun --mpi=pmix -n 64 pw.x -nk 4 -in s3_ref_aMn_noU.in > s3_ref_aMn_noU.out 2>&1
if grep -q "JOB DONE" s3_ref_aMn_noU.out; then
    touch s3_ref_aMn_noU.done
    echo "[ok] s3_ref_aMn_noU  $(grep '!    total energy' s3_ref_aMn_noU.out | tail -1)"
else
    echo "[FAIL] s3_ref_aMn_noU -- died BEFORE or AFTER '!    total energy'?"
    grep -m2 -A3 'Error in routine' s3_ref_aMn_noU.out
fi
rm -rf ./tmp_s3_ref_aMn_noU
