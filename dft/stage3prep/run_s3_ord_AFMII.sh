#!/bin/bash
#SBATCH --job-name=s3_ord_AFMII
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --output=slurm_s3_ord_AFMII_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage3prep

if [ -f "s3_ord_AFMII.done" ]; then echo "[skip] s3_ord_AFMII already complete"; exit 0; fi
mkdir -p ./tmp_s3_ord_AFMII
srun --mpi=pmix -n 64 pw.x -nk 4 -in s3_ord_AFMII.in > s3_ord_AFMII.out 2>&1
if grep -q "JOB DONE" s3_ord_AFMII.out; then
    touch s3_ord_AFMII.done
    echo "[ok] s3_ord_AFMII  $(grep '!    total energy' s3_ord_AFMII.out | tail -1)"
else
    echo "[FAIL] s3_ord_AFMII -- died BEFORE or AFTER '!    total energy'?"
    grep -m2 -A3 'Error in routine' s3_ord_AFMII.out
fi
rm -rf ./tmp_s3_ord_AFMII
