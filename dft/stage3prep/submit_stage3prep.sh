#!/bin/bash
#SBATCH --job-name=MnSe-s3prep
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage3prep

for tag in s3_rs_AFMII_prim s3_ref_Se s3_ref_bccMn_U s3_ref_bccMn_noU s3_ord_FM s3_ord_AFMI s3_ord_AFMII s3_ref_aMn_U s3_ref_aMn_noU; do
    if [ -f "${tag}.done" ]; then echo "[skip] $tag"; continue; fi
    echo "=== $tag  $(date) ==="
    mkdir -p ./tmp_${tag}
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}.in > ${tag}.out 2>&1
    if grep -q "JOB DONE" ${tag}.out; then
        touch ${tag}.done
        echo "[ok] $tag  $(grep '!    total energy' ${tag}.out | tail -1)"
    else
        echo "[FAIL] $tag -- died BEFORE or AFTER '!    total energy'?"
        grep -m2 -A3 'Error in routine' ${tag}.out
    fi
    rm -rf ./tmp_${tag}
done
echo "ALL DONE $(date)"
