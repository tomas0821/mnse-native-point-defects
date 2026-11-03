#!/bin/bash
#SBATCH --job-name=MnSe-probe3
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --time=00:40:00
#SBATCH --output=slurm_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/stage1
mkdir -p ./tmp_uv
srun --mpi=pmix -n 32 pw.x -in probe_UV.in > probe_UV.out 2>&1
echo "=== U+V result ==="
grep -iE "Error in routine|not implemented|not supported" probe_UV.out | head -4
grep -iE "Hubbard energy|inter-site V|Forces acting on atoms" probe_UV.out | head -5
grep -q "JOB DONE" probe_UV.out && echo "U+V relax: JOB DONE -- FORCES WORK" || echo "U+V relax: FAILED"
rm -rf ./tmp_uv
