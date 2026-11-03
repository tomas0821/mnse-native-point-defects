#!/bin/bash
#SBATCH --job-name=MnSe-anti
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=16:00:00
#SBATCH --output=slurm_anti_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/anti
for M in 0p0 1p0 2p0 3p0 4p0 5p0; do
    t=anti_M${M}
    [ -f ${t}.done ] && { echo "[skip] $t"; continue; }
    mkdir -p ./tmp_${t}
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${t}.in > ${t}.out 2>&1
    if grep -q "JOB DONE" ${t}.out && ! grep -q "convergence NOT achieved" ${t}.out; then
        touch ${t}.done
        echo "[ok] $t  iters=$(grep -c 'iteration #' ${t}.out)  $(grep '^!    total energy' ${t}.out|tail -1)"
    else
        echo "[FAIL] $t  iters=$(grep -c 'iteration #' ${t}.out)"
    fi
    rm -rf ./tmp_${t}
done
echo "ALL DONE $(date)"
