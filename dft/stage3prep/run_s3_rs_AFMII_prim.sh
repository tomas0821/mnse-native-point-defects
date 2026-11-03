#!/bin/bash
#SBATCH --job-name=s3_rs_AFMII_prim
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00
#SBATCH --output=slurm_s3_rs_AFMII_prim_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage3prep

if [ -f "s3_rs_AFMII_prim.done" ]; then echo "[skip] s3_rs_AFMII_prim already complete"; exit 0; fi
mkdir -p ./tmp_s3_rs_AFMII_prim
srun --mpi=pmix -n 16 pw.x -nk 8 -in s3_rs_AFMII_prim.in > s3_rs_AFMII_prim.out 2>&1
if grep -q "JOB DONE" s3_rs_AFMII_prim.out; then
    touch s3_rs_AFMII_prim.done
    echo "[ok] s3_rs_AFMII_prim  $(grep '!    total energy' s3_rs_AFMII_prim.out | tail -1)"
else
    echo "[FAIL] s3_rs_AFMII_prim -- died BEFORE or AFTER '!    total energy'?"
    grep -m2 -A3 'Error in routine' s3_rs_AFMII_prim.out
fi
rm -rf ./tmp_s3_rs_AFMII_prim
