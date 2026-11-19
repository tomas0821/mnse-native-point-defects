#!/bin/bash
#SBATCH --job-name=s5-defscf
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --output=slurm_s5_defscf_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage5

# One job, four SCFs, densities KEPT. Needed for the Stage 5 bands step AND for
# the spin-resolved PDOS -- i.e. not contingent on the pristine go/no-go.
for d in V_Se V_Mn Se_Mn Mn_Se; do
    t=wzsc_${d}
    [ -f ${t}.scfdone ] && { echo "[skip] $t"; continue; }
    mkdir -p ./tmp_${t}
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${t}_scf.in > ${t}_scf.out 2>&1
    if grep -q "JOB DONE" ${t}_scf.out && ! grep -q "convergence NOT achieved" ${t}_scf.out; then
        touch ${t}.scfdone
        echo "[ok] $t E=$(grep '^!    total energy' ${t}_scf.out|tail -1|awk '{print $5}') tot_m=$(grep 'total magnetization' ${t}_scf.out|tail -1|awk '{print $4}')"
    else
        echo "[FAIL] $t"
    fi
    # density deliberately NOT deleted
done
echo "S5 DEFSCF DONE $(date)"
