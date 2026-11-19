#!/bin/bash
# Relax one structure in legs of under 24 h. A leg that hits max_seconds resubmits itself from its last geometry
# (a fresh relax, as in prod/restart_relax.sh: it needs only the .out file, at the cost of the BFGS Hessian). At most 4 legs.
# Usage: sbatch --job-name=sr-NAME run_relax.sh NAME [LEG]
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --time=24:00:00
#SBATCH --output=slurm_%x_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
export OMP_NUM_THREADS=1
cd "$SLURM_SUBMIT_DIR"
t=$1; leg=${2:-1}
if [ "$leg" -eq 1 ]; then
  cp "$t.in" "${t}_L1.in"
else
  prev="${t}_L$((leg-1))"
  sed '/^ATOMIC_POSITIONS/,$d' "$prev.in" > "${t}_L$leg.in"
  awk '/^ATOMIC_POSITIONS/{delete L; n=0; on=1; L[n++]=$0; next}
       on&&/^ *(Mn1|Mn2|Se) +-?[0-9.]+ +-?[0-9.]+ +-?[0-9.]+/{L[n++]=$0; next}
       on{on=0}
       END{for(i=0;i<n;i++) print L[i]}' "$prev.out" >> "${t}_L$leg.in"
  sed -n '/^K_POINTS/,$p' "$t.in" >> "${t}_L$leg.in"
fi
srun --mpi=pmix -n 64 pw.x -nk 4 -in "${t}_L$leg.in" > "${t}_L$leg.out" 2>&1
o="${t}_L$leg.out"
echo "leg $leg: $(grep -c '^!' "$o") steps, $(grep '^!' "$o" | tail -1), $(grep 'Total force' "$o" | tail -1)"
if grep -q "End of BFGS Geometry Optimization" "$o"; then
  echo "converged"; touch "$t.done"; rm -rf "tmp_$t"
elif grep -q "convergence NOT achieved\|Error in routine" "$o"; then
  echo "FAILED, not resubmitting"; grep -E "Error in routine|convergence NOT" "$o" | head -3
elif [ "$leg" -lt 4 ]; then
  rm -rf "tmp_$t"; sbatch --job-name="sr-$t" run_relax.sh "$t" $((leg+1))
else
  echo "stopped after 4 legs, not converged"
fi
