# Recreating the Stage 5 wavefunctions

The wavefunction and density directories on cluster scratch were deleted on 2026-11-06 to free about 2.6 TB (user-approved):
`stage5/tmp_bands_{V_Mn,V_Se,Se_Mn,Mn_Se}` and `stage5/tmp_wzsc`, `stage5/tmp_wzsc_{V_Mn,V_Se,Se_Mn,Mn_Se}`.

**Nothing in the paper needs them.** Figure 4, Table 3 and Table S14 come from the small files that were kept: `w_*.npz` (spectral weights),
`wzsc_*_bands.out` (eigenvalues), and the `s_*.npz` and `spread_*.npz` derived from them. Recreate the wavefunctions only to unfold in a different
way, for example to change the sublattice projection or look at other k-points, or to compute a PDOS.

## What each piece costs

Measured on `mpi`, 64 cores, `pw.x -nk 4`. All jobs fit the cluster's 24 h limit.

| Step | Script | Original job | Wall time | Disk |
|---|---|---|---|---|
| Pristine SCF | `run_stage5_pristine.sh` | 121046 | 3 h 45 m | with the next row |
| Pristine bands | same script (or `run_s5_bands.sh`) | 121354 | 8 h 44 m | 487 GB in `tmp_wzsc` |
| Four defect SCFs, one job | `run_s5_defscf.sh` | 121375 | 15 h 44 m (3 h 39 m to 4 h 28 m each) | 41-42 GB each in `tmp_wzsc_<d>` |
| V_Mn bands | `run_bands_V_Mn.sh` | 140272 | 9 h 11 m | 485 GB in `tmp_bands_V_Mn` |
| V_Se bands | `run_bands_V_Se.sh` | 140273 | 11 h 17 m | 487 GB |
| Se_Mn bands | `run_bands_Se_Mn.sh` | 140274 | 9 h 9 m | 485 GB |
| Mn_Se bands | `run_bands_Mn_Se.sh` | 140275 | 11 h 39 m | 489 GB |
| Unfolding, each cell | `run_unfold.sh` | 140276, 140321-5, 140344 | under 1 min (`shared`, 16 cores) | 389 KB |

Total: about 69 h of 64-core jobs and about 2.6 TB, run in parallel over about a day. Check `df -h /scratch/global` first.

For one defect only: its SCF (about 4 h, run it by hand or cut the loop in `run_s5_defscf.sh`) and its bands (about 10 h).

## Files

On the cluster they sit together in `/scratch/global/tomas.rojas_s/mnse/stage5/`, which still holds every input, script and log. If scratch has been wiped, restore them
from the project's `dft/stage5/`: the top level holds `wzsc_scf.in`, `run_stage5_pristine.sh`, `run_s5_bands.sh`, `run_s5_defscf.sh`, `run_unfold.sh` and `unfold_weights.py`;
`defbands/` holds `wzsc_bands.in`, `wzsc_<d>_scf.in`, `wzsc_<d>_bands.in`, `run_bands_<d>.sh` and `gen_defbands.py`. All must go into one flat folder on the cluster, because the scripts
`cd` there and read inputs by bare name. Checksums matched the local copies on 2026-11-06.

Environment: `module load intel/2024.2.1 qe/7.5-cpu`, pseudopotentials in `/home/tomas.rojas_s/pseudo` (`Mn.pbe-spn-kjpaw_psl.0.3.1.UPF`, `Se.pbe-dn-kjpaw_psl.1.0.0.UPF`), and for unfolding
`module load python/3.11.9-gcc-12.2.0-kywavjj py-numpy/2.1.1-gcc-12.2.0-epyzfd2`. **Slurm only, never the login node.** If the VPN is down, only the user can restart it.

## Steps

```bash
cd /scratch/global/tomas.rojas_s/mnse/stage5

# 1. Pristine: SCF and bands in one job, wavefunctions left in tmp_wzsc
sbatch run_stage5_pristine.sh

# 2. Defect SCFs. run_s5_defscf.sh SKIPS any defect with a .scfdone marker, and those markers still exist, so remove them first
rm -f wzsc_V_Mn.scfdone wzsc_V_Se.scfdone wzsc_Se_Mn.scfdone wzsc_Mn_Se.scfdone
sbatch run_s5_defscf.sh

# 3. Defect bands, after step 2. Each copies the SCF save into a fresh tmp_bands_<d> without its wavefunctions
for d in V_Mn V_Se Se_Mn Mn_Se; do sbatch --dependency=afterok:<SCF job id> run_bands_$d.sh; done

# 4. Unfold, after steps 1 and 3
sbatch run_unfold.sh tmp_wzsc/wzsc.save w_pristine.npz
for d in V_Mn V_Se Se_Mn Mn_Se; do sbatch run_unfold.sh tmp_bands_$d/wzsc_$d.save w_$d.npz; done
```

Step 1 overwrites `wzsc_scf.out` and `wzsc_bands.out`, and steps 2-4 overwrite their `.out` and `w_*.npz` files. Move the old ones aside if you want to compare. The bands
scripts stop with `[FAIL] no density` if their SCF hasn't finished, hence the dependency in step 3. Unfolding has to wait for steps 1 and 3 in the same way.

Then, on the local machine, after copying `w_*.npz` and `wzsc_*_bands.out` to `dft/stage5/`:

```bash
python3 peak_track.py w_pristine.npz wzsc_bands.out --dump s_pristine.npz
for d in V_Mn V_Se Se_Mn Mn_Se; do python3 peak_track.py w_$d.npz wzsc_${d}_bands.out --dump s_$d.npz; done
python3 make_stage5_table.py      # Table 3, also writes s_relaxed_u.npz and spread_*.npz
python3 make_stage5_figure.py     # Figure 4
python3 estimator_compare.py      # Table S14 rows
```

## What the geometries are

- `wzsc_scf.in`: pristine 72-atom cell at the ideal internal parameter u = 0.375. This is deliberate: it is the cell the unfolding was validated on, against the ideal-u primitive cell.
- `wzsc_<d>_scf.in`: relaxed defect geometries from `dft/prod/`, at 120/960 Ry and 3x3x3 k. Their energies are the production energies in `formation_energies.py`.
- Mn_Se is constrained at `tot_magnetization = 3`, and that carries into its bands input.
- `wzsc_*_bands.in`: 41 k-points along L1-Gamma-L2, written in supercell crystal coordinates so that each supercell K equals a primitive k (L1 is 1.5 0 1.0 here, primitive 0.5 0 0.5). 590 bands, `verbosity = 'high'`, which must stay in `&CONTROL`: in the wrong namelist it made the first pristine bands run fail.
- `wzsc_scf.in` has `tstress = .true.`, which ran fine here. The similar `fs128` run crashed in `stres_hub.f90` after its SCF, so if that happens, set it to `.false.`: the energy and density don't need it.

## Check that it worked

| Cell | SCF energy (Ry) | M | Peak odd part at t = 13, sigma = 0.03 eV (meV) | `w` sha256 start |
|---|---|---|---|---|
| Pristine | -24780.54010775 | 0.00 | 239.0 | c8289aee817b |
| V_Mn | -24493.87426264 | -5.25 | 225.6 | 8032c3e47aa5 |
| V_Se | -24378.46008823 | 0.00 | 234.8 | 56288a99d666 |
| Se_Mn | -24895.70752742 | -5.00 | 230.5 | e96f84f81ef5 |
| Mn_Se | -24664.70559843 | 3.00 | 237.4 | 573f65fd14fa |

- Energies should match to about 1e-6 Ry. The checksums may not match on a different build or core count, so compare numbers, not bits.
- Each bands log should say `k-blocks=82 wfcup=41`.
- Each unfolding log prints the total weight at the first k-point. Pristine must give `sum_w(up,k1)=33.000 sum_w(dn,k1)=33.000`; V_Mn gave 33.426 and 32.791.
- On the pristine cell, the tracked curve must reproduce the ideal-u primitive curve (`uideal_bands.out`) to about 1.5 meV. Table S10 lists the rest.

Original logs are in `recreate_logs/`.

## If you want a PDOS instead

The defect SCF step is enough; the bands and unfolding steps are not needed. Run the SCF of step 2, then an `nscf` on a denser grid and `projwfc.x`, in the same `outdir`. The nscf input hasn't been written. Budget about 4 h per defect for the SCF, plus the nscf.
