#!/bin/bash
#SBATCH --job-name=d_rs_V_Se_s120
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=08:00:00
#SBATCH --output=slurm_d_rs_V_Se_s120_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/prod
awk '/Begin final coordinates/,/End final coordinates/' ../pilot/p_VSe_070.out \
  | sed -n '/CELL_PARAMETERS/,$p' | grep -v "End final coordinates" > geom_d_rs_V_Se.txt
[ -s geom_d_rs_V_Se.txt ] || { echo "[FAIL] no relaxed geometry"; exit 1; }
cat d_rs_V_Se_s120.head geom_d_rs_V_Se.txt d_rs_V_Se_kpts.txt > d_rs_V_Se_s120.in
mkdir -p ./tmp_d_rs_V_Se_s120
srun --mpi=pmix -n 64 pw.x -nk 4 -in d_rs_V_Se_s120.in > d_rs_V_Se_s120.out 2>&1
grep -q "JOB DONE" d_rs_V_Se_s120.out && ! grep -q "convergence NOT achieved" d_rs_V_Se_s120.out \
  && { touch d_rs_V_Se.done; echo "[ok] $(grep '^!    total energy' d_rs_V_Se_s120.out|tail -1)"; } \
  || echo "[FAIL] d_rs_V_Se_s120"
rm -rf ./tmp_d_rs_V_Se_s120
