#!/bin/bash
#SBATCH --job-name=MnSe-stage1
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu

cd /scratch/global/tomas.rojas_s/mnse/stage1

for tag in convA_rs_ecut40 convA_rs_ecut50 convA_rs_ecut60 convA_rs_ecut70 convA_rs_ecut80 convA_rs_ecut90 convA_wz_ecut40 convA_wz_ecut50 convA_wz_ecut60 convA_wz_ecut70 convA_wz_ecut80 convA_wz_ecut90 convA_rs_k4 convA_rs_k6 convA_rs_k8 convA_rs_k10 convA_wz_k6 convA_wz_k8 convA_wz_k10 convA_wz_k12 ordB_rs_FM ordB_rs_AFMI ordB_rs_AFMII ordB_wz_FM ordB_wz_AFM relC_wz_AFM relC_rs_AFMI_quick relC_rs_AFMII; do
    if [ -f "${tag}.done" ]; then
        echo "[skip] $tag already done"
        continue
    fi
    echo "=== $tag  $(date) ==="
    # mkdir BEFORE the run. The loop deletes ./tmp at the end of each iteration,
    # and 64 MPI ranks then race to recreate it -- some lose and pw.x dies with
    # "unable to create directory ./tmp/<prefix>.save/" AFTER the SCF has already
    # converged. Cost us a 4 h 64-atom job on 2026-08-07. Creating it up front,
    # serially, removes the race.
    mkdir -p ./tmp
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
