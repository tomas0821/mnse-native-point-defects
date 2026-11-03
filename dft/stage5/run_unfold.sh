#!/bin/bash
#SBATCH --job-name=s5-unfold
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=6:00:00
#SBATCH --output=slurm_unfold_%j.out
# usage: sbatch run_unfold.sh <save_dir> <out.npz>
module purge
module load python/3.11.9-gcc-12.2.0-kywavjj py-numpy/2.1.1-gcc-12.2.0-epyzfd2
cd /scratch/global/tomas.rojas_s/mnse/stage5
python3 -c "import numpy; print('numpy', numpy.__version__)" || { echo "[FAIL] numpy"; exit 1; }
echo "start $(date)"
export OMP_NUM_THREADS=1
srun --exact -n 1 -c $SLURM_CPUS_PER_TASK python3 -u unfold_weights.py "$1" "$2" --M 3 3 2 --nproc $SLURM_CPUS_PER_TASK && echo "[ok] $2" || echo "[FAIL] unfold"
echo "UNFOLD DONE $(date)"
