#!/bin/bash
#SBATCH --job-name=MnSe-restart
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --output=slurm_restart_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/prod
t=d_wz_Se_Mn
bash restart_relax.sh $t
grep -cE '^ *(Mn1|Mn2|Se) ' geom_rs_${t}.txt
mkdir -p ./tmp_${t}_rs
srun --mpi=pmix -n 64 pw.x -nk 4 -in ${t}_rs.in > ${t}_rs.out 2>&1
if grep -q "End of BFGS Geometry Optimization" ${t}_rs.out; then
  echo "[RELAXED] $t  $(grep -c '^!    total energy' ${t}_rs.out) steps  $(grep '^!    total energy' ${t}_rs.out|tail -1)"
else
  echo "[INCOMPLETE] $t  $(grep -c '^!    total energy' ${t}_rs.out) steps  force=$(grep 'Total force' ${t}_rs.out|tail -1|awk '{print $4}')"
fi
