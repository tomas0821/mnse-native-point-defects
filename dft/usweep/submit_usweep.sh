#!/bin/bash
#SBATCH --job-name=MnSe-usweep
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu

cd /scratch/global/tomas.rojas_s/mnse/usweep

for tag in uA_wz_U0p0 uA_wz_U2p0 uA_wz_U3p0 uA_wz_U4p0 uA_wz_U5p1 uA_wz_U6p0 uB_wz_L_U6p0_J0p9 uA_rs_U0p0 uA_rs_U2p0 uA_rs_U3p0 uA_rs_U4p0 uA_rs_U5p1 uA_rs_U6p0 uB_rs_L_U6p0_J0p9; do
    if [ -f "${tag}.done" ]; then echo "[skip] $tag"; continue; fi
    echo "=== $tag  $(date) ==="
    mkdir -p ./tmp          # serial mkdir: avoids the 64-rank race that killed
                            # ordB_rs_FM_c120 after its SCF had already converged
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}.in > ${tag}.out 2>&1
    if grep -q "JOB DONE" ${tag}.out; then
        touch ${tag}.done
        echo "[ok] $tag  $(grep '!    total energy' ${tag}.out | tail -1)"
    else
        echo "[FAIL] $tag -- check whether it died BEFORE or AFTER '!    total energy'"
    fi
    rm -rf ./tmp
done
echo "ALL DONE $(date)"
