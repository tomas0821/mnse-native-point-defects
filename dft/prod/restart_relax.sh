#!/bin/bash
# Restart a wall-clock-killed relax from its LAST written geometry.
# A fresh relax (not restart_mode='restart') is used deliberately: it only needs
# the .out file, which always survives, whereas restart_mode depends on the
# scratch tmp_ dir that the job scripts may have removed. The cost is the lost
# BFGS Hessian -- a few extra ionic steps -- which is cheap insurance.
t=$1
cd /scratch/global/tomas.rojas_s/mnse/prod
sed -n '/^CELL_PARAMETERS/,/^ATOMIC_POSITIONS/p' ${t}_r70.in | head -4 > geom_rs_${t}.txt
awk '/^ATOMIC_POSITIONS/{blk=NR} {L[NR]=$0} END{
   for(i=blk;i<=NR;i++){ if(L[i]~/^K_POINTS/) break; if(L[i]~/^ATOMIC_POSITIONS/||L[i]~/^ *(Mn1|Mn2|Se) /) print L[i] }
}' ${t}_r70.out >> geom_rs_${t}.txt
n=$(grep -cE '^ *(Mn1|Mn2|Se) ' geom_rs_${t}.txt)
echo "  $t: restarting from $n atoms, last force $(grep 'Total force' ${t}_r70.out|tail -1|awk '{print $4}')"
sed '/^CELL_PARAMETERS/,$d' ${t}_r70.in | sed "s/prefix      = '${t}_r70'/prefix      = '${t}_rs'/; s|outdir      = './tmp_${t}_r70'|outdir      = './tmp_${t}_rs'|" > ${t}_rs.head
cat ${t}_rs.head geom_rs_${t}.txt ${t}_kpts.txt > ${t}_rs.in
