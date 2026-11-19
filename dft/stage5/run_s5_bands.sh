#!/bin/bash
#SBATCH --job-name=s5-bands
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --output=slurm_s5_bands_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage5

# The SCF density from job 121046 survived in ./tmp_wzsc -- reuse it, do NOT redo
# the 13-iteration SCF. Only the bands step failed (verbosity was in the wrong
# namelist).
[ -f tmp_wzsc/wzsc.save/charge-density.dat ] || { echo "[FAIL] no charge density"; exit 1; }
srun --mpi=pmix -n 64 pw.x -nk 4 -in wzsc_bands.in > wzsc_bands.out 2>&1
if grep -q "JOB DONE" wzsc_bands.out && ! grep -qi "Error in routine" wzsc_bands.out; then
    echo "[ok] bands  k-blocks=$(grep -c '^          k =' wzsc_bands.out)"
else
    echo "[FAIL] bands"; grep -A3 "Error in routine" wzsc_bands.out | head -6
fi
echo "S5 BANDS DONE $(date)"
