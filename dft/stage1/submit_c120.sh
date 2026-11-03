#!/bin/bash
#SBATCH --job-name=MnSe-c120
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --output=slurm_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu

cd /scratch/global/tomas.rojas_s/mnse/stage1

for tag in ordB_rs_AFMII_c120 ordB_rs_AFMI_c120 ordB_rs_FM_c120 ordB_wz_AFM_c120 ordB_wz_FM_c120 relC_rs_AFMI_quick_c120 relC_wz_AFM_c120 ; do
    if [ -f "${tag}.done" ]; then
        echo "[skip] $tag already done"
        continue
    fi
    echo "=== $tag  $(date) ==="
    # This cluster PROHIBITS mpirun/mpiexec (SLURMMASTER rejects the job at
    # submit time). srun --mpi=pmix is mandatory. The CrN gate1 script used
    # mpirun and predates this policy -- do not copy it verbatim again.
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}.in > ${tag}.out 2>&1
    if grep -q "JOB DONE" ${tag}.out; then
        touch ${tag}.done
        echo "[ok] $tag  $(grep '!    total energy' ${tag}.out | tail -1)"
    else
        echo "[FAIL] $tag -- see ${tag}.out"
    fi
    rm -rf ./tmp
done
echo "ALL DONE $(date)"
