#!/bin/bash
#SBATCH --job-name=s5b-V_Mn
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --output=slurm_s5b_V_Mn_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage5
# Stage 5 defect bands (L1-Gamma-L2, 41 k, nbnd 590) on the job-121375 SCF density.
# Work in a FRESH outdir holding a copy of the SCF save minus its wavefunctions, so the
# SCF save in tmp_wzsc_V_Mn is never overwritten. Wavefunctions are KEPT for unfolding.
P=wzsc_V_Mn
[ -f tmp_${P}/${P}.save/charge-density.dat ] || { echo "[FAIL] no density"; exit 1; }
mkdir -p tmp_bands_V_Mn/${P}.save
find tmp_${P}/${P}.save -maxdepth 1 -type f ! -name 'wfc*.dat' -exec cp {} tmp_bands_V_Mn/${P}.save/ \;
srun --mpi=pmix -n 64 pw.x -nk 4 -in ${P}_bands.in > ${P}_bands.out 2>&1
if grep -q "JOB DONE" ${P}_bands.out && ! grep -qi "Error in routine" ${P}_bands.out; then
    echo "[ok] ${P} bands k-blocks=$(grep -c '^          k =' ${P}_bands.out) wfcup=$(ls tmp_bands_V_Mn/${P}.save | grep -c wfcup)"
else
    echo "[FAIL] ${P} bands"; grep -A3 "Error in routine" ${P}_bands.out | head -6
fi
echo "S5B V_Mn DONE $(date)"
