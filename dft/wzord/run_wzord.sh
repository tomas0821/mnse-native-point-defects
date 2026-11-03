#!/bin/bash
#SBATCH --job-name=MnSe-wzord
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --time=01:00:00
#SBATCH --output=slurm_wzord_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/wzord
for tag in wz376_AFM wz376_FM; do
    mkdir -p ./tmp_${tag}
    srun --mpi=pmix -n 16 pw.x -nk 8 -in ${tag}.in > ${tag}.out 2>&1
    grep -q "JOB DONE" ${tag}.out && echo "[ok] $tag  $(grep '!    total energy' ${tag}.out|tail -1)" \
                                  || echo "[FAIL] $tag"
    rm -rf ./tmp_${tag}
done
echo ALL DONE
