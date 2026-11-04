#!/bin/bash
# Constrained-magnetization single points at the RELAXED defect geometry, 70/560 Ry.
# Usage: sbatch run_scan.sh TAG M1 M2 ...   e.g. sbatch run_scan.sh d_rs_Se_Mn 5 1 7 3
#SBATCH --job-name=spinscan
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=16:00:00
#SBATCH --output=slurm_%x_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
export OMP_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
t=$1; shift
for M in "$@"; do
  n=${t}_M${M}
  sed -e "s/ecutwfc *= *120.0/ecutwfc     = 70.0/" -e "s/ecutrho *= *960.0/ecutrho     = 560.0/" \
      -e "s/^ *prefix *=.*/    prefix      = '$n'/" -e "s#^ *outdir *=.*#    outdir      = './tmp_$n'#" \
      -e "s/tstress *= *.true./tstress = .false./" -e "s/tprnfor *= *.true./tprnfor = .false./" \
      -e "s/^ *nspin *= *2/    nspin       = 2\n    tot_magnetization = $M/" \
      -e "s/mixing_beta *=.*/mixing_beta     = 0.1\n    mixing_mode     = 'local-TF'/" \
      -e "s/electron_maxstep *=.*/electron_maxstep = 500/" ../prod/${t}_s120.in > $n.in
  srun --mpi=pmix -n 64 pw.x -nk 4 -in $n.in > $n.out 2>&1
  echo "$n $(grep -c 'iteration #' $n.out) iters: $(grep '^!' $n.out | tail -1) | conv-failed $(grep -c 'convergence NOT achieved' $n.out)"
  rm -rf tmp_$n
done
