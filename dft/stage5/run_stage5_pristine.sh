#!/bin/bash
#SBATCH --job-name=s5-wzsc
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --output=slurm_s5_pristine_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage5

mkdir -p ./tmp_wzsc
srun --mpi=pmix -n 64 pw.x -nk 4 -in wzsc_scf.in > wzsc_scf.out 2>&1
if ! grep -q "JOB DONE" wzsc_scf.out || grep -q "convergence NOT achieved" wzsc_scf.out; then
    echo "[FAIL] scf"; exit 1
fi
echo "[ok] scf $(grep '^!    total energy' wzsc_scf.out | tail -1)"

# bands reads the density from ./tmp_wzsc -- do NOT delete it between the two steps
srun --mpi=pmix -n 64 pw.x -nk 4 -in wzsc_bands.in > wzsc_bands.out 2>&1
if grep -q "JOB DONE" wzsc_bands.out; then
    echo "[ok] bands  kpts=$(grep -c '^          k =' wzsc_bands.out)"
else
    echo "[FAIL] bands"
fi
echo "STAGE5 PRISTINE DONE $(date)"
