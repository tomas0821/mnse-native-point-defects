#!/bin/bash
#SBATCH --job-name=d_rs_Mn_Se
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --output=slurm_d_rs_Mn_Se_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/prod

# --- 1. relax at 70 Ry (forces shown to match 120 Ry to <0.1% in the V_Se pilot)
mkdir -p ./tmp_d_rs_Mn_Se_r70
srun --mpi=pmix -n 64 pw.x -nk 4 -in d_rs_Mn_Se_r70.in > d_rs_Mn_Se_r70.out 2>&1
if ! grep -q "JOB DONE" d_rs_Mn_Se_r70.out || grep -q "convergence NOT achieved" d_rs_Mn_Se_r70.out; then
    echo "[INCOMPLETE relax] d_rs_Mn_Se -- $(grep -c '^!    total energy' d_rs_Mn_Se_r70.out) ionic steps done"
    rm -rf ./tmp_d_rs_Mn_Se_r70; exit 0
fi
echo "[ok relax] d_rs_Mn_Se  $(grep -c '^!    total energy' d_rs_Mn_Se_r70.out) steps  $(grep '^!    total energy' d_rs_Mn_Se_r70.out|tail -1)"
rm -rf ./tmp_d_rs_Mn_Se_r70

# --- 2. single point at 120 Ry on the relaxed geometry (matches the mu references)
awk '/Begin final coordinates/,/End final coordinates/' d_rs_Mn_Se_r70.out \
  | sed -n '/CELL_PARAMETERS/,$p' | grep -v "End final coordinates" > geom_d_rs_Mn_Se.txt
if [ ! -s geom_d_rs_Mn_Se.txt ]; then echo "[FAIL] no final geometry for d_rs_Mn_Se"; exit 1; fi
cat d_rs_Mn_Se_s120.head geom_d_rs_Mn_Se.txt d_rs_Mn_Se_kpts.txt > d_rs_Mn_Se_s120.in
mkdir -p ./tmp_d_rs_Mn_Se_s120
srun --mpi=pmix -n 64 pw.x -nk 4 -in d_rs_Mn_Se_s120.in > d_rs_Mn_Se_s120.out 2>&1
grep -q "JOB DONE" d_rs_Mn_Se_s120.out && ! grep -q "convergence NOT achieved" d_rs_Mn_Se_s120.out && { touch d_rs_Mn_Se.done; echo "[ok 120Ry] d_rs_Mn_Se  $(grep '^!    total energy' d_rs_Mn_Se_s120.out|tail -1)"; } \
                                 || echo "[FAIL 120Ry] d_rs_Mn_Se"
rm -rf ./tmp_d_rs_Mn_Se_s120
