#!/bin/bash
#SBATCH --job-name=wzMnSe_M3rs
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --output=slurm_wzMnSe_M3rs_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/prod
# resume the constrained wurtzite antisite from its last written geometry
sed -n '/^CELL_PARAMETERS/,/^ATOMIC_POSITIONS/p' wzMnSe_M3.in | head -4 > geom_wzM3rs.txt
awk '/^ATOMIC_POSITIONS/{b=NR} {L[NR]=$0} END{
  for(i=b;i<=NR;i++){ if(L[i]~/^K_POINTS/) break
    if(L[i]~/^ATOMIC_POSITIONS/||L[i]~/^ *(Mn1|Mn2|Se) /) print L[i] } }' wzMnSe_M3.out >> geom_wzM3rs.txt
n=$(grep -cE '^ *(Mn1|Mn2|Se) ' geom_wzM3rs.txt)
echo "resuming from $n atoms, last force $(grep 'Total force' wzMnSe_M3.out|tail -1|awk '{print $4}')"
[ "$n" -eq 72 ] || { echo "[FAIL] expected 72 atoms, got $n"; exit 1; }
sed '/^CELL_PARAMETERS/,$d' wzMnSe_M3.in \
  | sed "s/prefix      = 'wzMnSe_M3'/prefix      = 'wzMnSe_M3rs'/; s|outdir      = './tmp_d_wz_Mn_Se_r70'|outdir      = './tmp_wzMnSe_M3rs'|" > wzM3rs.head
cat wzM3rs.head geom_wzM3rs.txt wzMnSe_kpts.txt > wzMnSe_M3rs.in
mkdir -p ./tmp_wzMnSe_M3rs
srun --mpi=pmix -n 64 pw.x -nk 4 -in wzMnSe_M3rs.in > wzMnSe_M3rs.out 2>&1
if grep -q "End of BFGS Geometry Optimization" wzMnSe_M3rs.out; then
  echo "[RELAXED] wzMnSe_M3rs $(grep -c '^!    total energy' wzMnSe_M3rs.out) steps $(grep '^!    total energy' wzMnSe_M3rs.out|tail -1)"
else
  echo "[INCOMPLETE] $(grep -c '^!    total energy' wzMnSe_M3rs.out) steps force=$(grep 'Total force' wzMnSe_M3rs.out|tail -1|awk '{print $4}')"
fi
