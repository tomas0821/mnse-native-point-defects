#!/bin/bash
#SBATCH --job-name=MnSe-s120
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=12:00:00
#SBATCH --output=slurm_s120batch_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/prod

# srun MUST be in this top-level script -- SLURMMASTER rejects the job if the
# only srun lives in a helper script it cannot see.
for t in d_rs_V_Se d_rs_Se_Mn d_rs_V_Mn d_rs_Mn_Se d_wz_V_Se d_wz_V_Mn d_wz_Mn_Se d_wz_Se_Mn; do
    [ -f ${t}.done ] && { echo "[skip] $t already done"; continue; }
    # a restarted relax wrote <tag>_rs.out; prefer it over the original
    src=${t}_r70.out; [ -f ${t}_rs.out ] && grep -q "End of BFGS" ${t}_rs.out && src=${t}_rs.out
    [ -f "$src" ] || { echo "[skip] $t: no relax output yet"; continue; }
    grep -q "End of BFGS Geometry Optimization" $src || { echo "[skip] $t: relax not converged yet"; continue; }

    # cell from the INPUT (fixed-cell relax never reprints it); positions from the output
    sed -n '/^CELL_PARAMETERS/,/^ATOMIC_POSITIONS/p' ${t}_r70.in | head -4 > geom_${t}.txt
    awk '/Begin final coordinates/,/End final coordinates/' $src \
      | sed -n '/^ATOMIC_POSITIONS/,$p' | grep -v "End final coordinates" >> geom_${t}.txt
    n=$(grep -cE "^\s*(Mn1|Mn2|Se)\s" geom_${t}.txt)
    echo "=== $t : $n atoms extracted ==="
    [ "$n" -gt 0 ] || { echo "[FAIL] $t extraction empty"; continue; }

    cat ${t}_s120.head geom_${t}.txt ${t}_kpts.txt > ${t}_s120.in
    mkdir -p ./tmp_${t}_s120
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${t}_s120.in > ${t}_s120.out 2>&1
    if grep -q "JOB DONE" ${t}_s120.out && ! grep -q "convergence NOT achieved" ${t}_s120.out; then
        touch ${t}.done
        echo "[ok] $t  $(grep '^!    total energy' ${t}_s120.out | tail -1)"
    else
        echo "[FAIL] $t 120Ry single point"
    fi
    rm -rf ./tmp_${t}_s120
done
echo "ALL DONE $(date)"
