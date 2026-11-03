#!/bin/bash
#SBATCH --job-name=d_wz_V_Mn_rs
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --output=slurm_d_wz_V_Mn_rs_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/prod
bash restart_relax.sh d_wz_V_Mn
mkdir -p ./tmp_d_wz_V_Mn_rs
srun --mpi=pmix -n 64 pw.x -nk 4 -in d_wz_V_Mn_rs.in > d_wz_V_Mn_rs.out 2>&1
if grep -q "End of BFGS Geometry Optimization" d_wz_V_Mn_rs.out; then
  echo "[RELAXED] d_wz_V_Mn $(grep -c '^!    total energy' d_wz_V_Mn_rs.out) steps $(grep '^!    total energy' d_wz_V_Mn_rs.out|tail -1)"
else
  echo "[INCOMPLETE] d_wz_V_Mn $(grep -c '^!    total energy' d_wz_V_Mn_rs.out) steps force=$(grep 'Total force' d_wz_V_Mn_rs.out|tail -1|awk '{print $4}')"
fi
