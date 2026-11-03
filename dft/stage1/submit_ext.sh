#!/bin/bash
#SBATCH --job-name=MnSe-ext
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=slurm_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu

cd /scratch/global/tomas.rojas_s/mnse/stage1

for tag in convA_rs_ecut100 convA_rs_ecut110 convA_rs_ecut120 convA_rs_ecut130 convA_rs_ecut140 convA_wz_ecut100 convA_wz_ecut110 convA_wz_ecut120 convA_wz_ecut130 convA_wz_ecut140 ; do
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
