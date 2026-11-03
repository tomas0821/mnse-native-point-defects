#!/bin/bash
#SBATCH --job-name=MnSe-probe2
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=01:00:00
#SBATCH --output=slurm_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage1

echo "############ PROBE 2: DFT+U+V forces ############"
mkdir -p ./tmp_uv
srun --mpi=pmix -n 32 pw.x -in probe_UV.in > probe_UV.out 2>&1
grep -iE "Hubbard parameters|inter-site|V\(|Error in routine|not implemented|Forces acting" probe_UV.out | head -10
grep -cE "^!    total energy" probe_UV.out | sed 's/^/scf energies: /'
grep -q "JOB DONE" probe_UV.out && echo "UV: JOB DONE (forces OK)" || echo "UV: FAILED"
rm -rf ./tmp_uv

echo "############ PROBE 1: hp.x linear response ############"
mkdir -p ./tmp_hp
srun --mpi=pmix -n 32 pw.x -in probe_hp_scf.in > probe_hp_scf.out 2>&1
grep -q "JOB DONE" probe_hp_scf.out && echo "scf OK" || { echo "scf FAILED"; tail -5 probe_hp_scf.out; }
srun --mpi=pmix -n 32 hp.x -in probe_hp.in > probe_hp.out 2>&1
grep -iE "Error in routine|not implemented|not supported|Hubbard U|chi0|converged" probe_hp.out | head -12
grep -q "JOB DONE" probe_hp.out && echo "hp.x: JOB DONE" || echo "hp.x: FAILED"
rm -rf ./tmp_hp
echo "ALL PROBES DONE"
