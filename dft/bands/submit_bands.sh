#!/bin/bash
#SBATCH --job-name=MnSe-bands
#SBATCH --partition=mpi
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=08:00:00
#SBATCH --output=slurm_%j.out

module purge
module load intel/2024.2.1 2>/dev/null
module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/bands

# splice the relaxed cell+positions out of a vc-relax output
extract_geom () {   # $1 = vc-relax .out
    awk '/Begin final coordinates/,/End final coordinates/' "$1" \
      | sed -n '/CELL_PARAMETERS/,$p' \
      | grep -v "End final coordinates"
}

for tag in b_wz b_rs; do
    echo "######## $tag  $(date) ########"
    mkdir -p ./tmp_${tag}

    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}_relax.in > ${tag}_relax.out 2>&1
    if ! grep -q "JOB DONE" ${tag}_relax.out; then echo "  [FAIL] relax"; continue; fi
    echo "  relax: $(grep '^!    total energy' ${tag}_relax.out | tail -1)"

    extract_geom ${tag}_relax.out > geom_${tag}.txt
    if [ ! -s geom_${tag}.txt ]; then echo "  [FAIL] no final geometry"; continue; fi

    cat ${tag}_scf.head   geom_${tag}.txt > ${tag}_scf.in
    printf 'K_POINTS automatic\n' >> ${tag}_scf.in
    grep -A1 "K_POINTS automatic" ${tag}_relax.in | tail -1 >> ${tag}_scf.in

    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}_scf.in > ${tag}_scf.out 2>&1
    if grep -q "JOB DONE" ${tag}_scf.out; then
        echo "  GAP:"; grep -iE "highest occupied|lowest unoccupied" ${tag}_scf.out | tail -2
    else
        echo "  [FAIL] scf"; grep -A3 'Error in routine' ${tag}_scf.out | head -6; continue
    fi

    cat ${tag}_bands.head geom_${tag}.txt ${tag}_kpath.txt > ${tag}_bands.in
    srun --mpi=pmix -n 64 pw.x -nk 4 -in ${tag}_bands.in > ${tag}_bands.out 2>&1
    grep -q "JOB DONE" ${tag}_bands.out && echo "  [ok] bands done" || echo "  [FAIL] bands"
done
echo "ALL DONE $(date)"
