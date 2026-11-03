#!/bin/bash
#SBATCH --job-name=MnSe-sym
#SBATCH --partition=shared
#SBATCH --account=c4139
#SBATCH --qos=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --time=03:00:00
#SBATCH --output=slurm_sym_%j.out
module purge; module load intel/2024.2.1 2>/dev/null; module load qe/7.5-cpu
cd /scratch/global/tomas.rojas_s/mnse/sym

echo "######## A. NULL CONTROL: alpha-MnSe AFM-II primitive, must give ZERO splitting ########"
mkdir -p ./tmp_null_rs
srun --mpi=pmix -n 16 pw.x -nk 4 -in null_rs_scf.in > null_rs_scf.out 2>&1
grep -q "JOB DONE" null_rs_scf.out && echo "  scf ok" || { echo "  scf FAILED"; grep -m1 -A3 'Error in routine' null_rs_scf.out; }
srun --mpi=pmix -n 16 pw.x -nk 4 -in null_rs_bands.in > null_rs_bands.out 2>&1
grep -q "JOB DONE" null_rs_bands.out && echo "  bands ok ($(grep -c 'bands (ev)' null_rs_bands.out) blocks)" || {
    echo "  bands FAILED -- retrying with cg diagonalisation"
    sed 's|    mixing_beta|    diagonalization = '"'"'cg'"'"'\n    mixing_beta|' null_rs_bands.in > null_rs_bands_cg.in
    srun --mpi=pmix -n 16 pw.x -nk 4 -in null_rs_bands_cg.in > null_rs_bands.out 2>&1
    grep -q "JOB DONE" null_rs_bands.out && echo "  bands ok with cg" || echo "  still FAILED"
}

echo "######## B. g-WAVE SCAN: w-MnSe azimuthal k-circle ########"
mkdir -p ./tmp_gw
# reuse the relaxed w-MnSe from the bands run
cp ../bands/b_wz_scf.in gw_scf.in 2>/dev/null
sed -i "s/prefix      = 'b_wz'/prefix      = 'gw'/; s|outdir      = './tmp_b_wz'|outdir      = './tmp_gw'|" gw_scf.in
srun --mpi=pmix -n 16 pw.x -nk 4 -in gw_scf.in > gw_scf.out 2>&1
grep -q "JOB DONE" gw_scf.out && echo "  scf ok" || { echo "  scf FAILED"; grep -m1 -A3 'Error in routine' gw_scf.out; }
sed -e "s/calculation = 'scf'/calculation = 'bands'/" -e "/&CONTROL/a\\    verbosity   = 'high'" gw_scf.in > gw_bands_head.in
sed '/K_POINTS/,$d' gw_bands_head.in > gw_bands.in
cat gwave_kpts.txt >> gw_bands.in
srun --mpi=pmix -n 16 pw.x -nk 4 -in gw_bands.in > gw_bands.out 2>&1
grep -q "JOB DONE" gw_bands.out && echo "  g-wave scan ok ($(grep -c 'bands (ev)' gw_bands.out) blocks)" || echo "  g-wave scan FAILED"
rm -rf ./tmp_null_rs ./tmp_gw
echo "ALL DONE"
