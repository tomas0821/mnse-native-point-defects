#!/bin/bash
#SBATCH --job-name=rsMnSe_M5_s120
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=10:00:00
#SBATCH --output=slurm_rsM5_s120_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/anti
sed -n '/^CELL_PARAMETERS/,/^ATOMIC_POSITIONS/p' rsMnSe_M5.in | head -4 > geom_rsM5.txt
awk '/Begin final coordinates/,/End final coordinates/' rsMnSe_M5.out \
  | sed -n '/^ATOMIC_POSITIONS/,$p' | grep -v "End final coordinates" >> geom_rsM5.txt
n=$(grep -cE '^ *(Mn1|Mn2|Se) ' geom_rsM5.txt); echo "geometry: $n atoms"
[ "$n" -eq 64 ] || { echo "[FAIL] expected 64, got $n"; exit 1; }
sed '/^CELL_PARAMETERS/,$d' rsMnSe_M5.in \
  | sed "s/prefix      = 'rsMnSe_M5'/prefix      = 'rsMnSe_M5_s120'/; s|outdir      = './tmp_rsMnSe_M5'|outdir      = './tmp_rsM5s'|
         s/calculation = 'relax'/calculation = 'scf'/; s/ecutwfc     = 70.0/ecutwfc     = 120.0/; s/ecutrho     = 560.0/ecutrho     = 960.0/" \
  | sed '/&IONS/,/^\//d' | sed '/forc_conv_thr\|etot_conv_thr\|nstep/d' > rsM5s.head
printf 'K_POINTS automatic\n  3 3 3 0 0 0\n' > rsM5_kpts.txt
cat rsM5s.head geom_rsM5.txt rsM5_kpts.txt > rsMnSe_M5_s120.in
mkdir -p ./tmp_rsM5s
srun --mpi=pmix -n 64 pw.x -nk 4 -in rsMnSe_M5_s120.in > rsMnSe_M5_s120.out 2>&1
grep -q "JOB DONE" rsMnSe_M5_s120.out && ! grep -q "convergence NOT" rsMnSe_M5_s120.out \
  && echo "[ok] $(grep '^!    total energy' rsMnSe_M5_s120.out|tail -1)" || echo "[FAIL]"
rm -rf ./tmp_rsM5s
