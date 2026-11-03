#!/bin/bash
# Build + run the 120 Ry single point from a finished fixed-cell relax.
# A fixed-cell `relax` prints "Begin final coordinates" WITHOUT CELL_PARAMETERS
# (the cell never changed), so the cell must come from the INPUT file and only
# the positions from the output. Extracting both from the output yields nothing.
t=$1
cd /scratch/global/tomas.rojas_s/mnse/prod
[ -f ${t}_r70.out ] || { echo "[skip] $t: no relax output"; exit 0; }
grep -q "End of BFGS Geometry Optimization" ${t}_r70.out || { echo "[skip] $t: relax not converged"; exit 0; }

# cell from the input (fixed during relax)
sed -n '/^CELL_PARAMETERS/,/^ATOMIC_POSITIONS/p' ${t}_r70.in | head -4 > geom_${t}.txt
# final positions from the output
awk '/Begin final coordinates/,/End final coordinates/' ${t}_r70.out \
  | sed -n '/^ATOMIC_POSITIONS/,$p' | grep -v "End final coordinates" >> geom_${t}.txt
n=$(grep -cE "^\s*(Mn1|Mn2|Se)\s" geom_${t}.txt)
echo "  $t: geometry has $n atoms"
[ "$n" -gt 0 ] || { echo "[FAIL] $t: extraction still empty"; exit 1; }
cat ${t}_s120.head geom_${t}.txt ${t}_kpts.txt > ${t}_s120.in
mkdir -p ./tmp_${t}_s120
srun --mpi=pmix -n 64 pw.x -nk 4 -in ${t}_s120.in > ${t}_s120.out 2>&1
if grep -q "JOB DONE" ${t}_s120.out && ! grep -q "convergence NOT achieved" ${t}_s120.out; then
    touch ${t}.done; echo "[ok] $t  $(grep '^!    total energy' ${t}_s120.out|tail -1)"
else
    echo "[FAIL] $t 120Ry single point"
fi
rm -rf ./tmp_${t}_s120
