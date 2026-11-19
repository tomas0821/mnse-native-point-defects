# Native point defects in rock-salt and wurtzite MnSe — calculation data

Data, inputs and analysis scripts for a first-principles study of the neutral native point
defects (V<sub>Mn</sub>, V<sub>Se</sub>, Se<sub>Mn</sub>, Mn<sub>Se</sub>) in the two MnSe
polymorphs: rock-salt α-MnSe (type-II antiferromagnet) and wurtzite w-MnSe.

**On the magnetic order of w-MnSe.** The defect cells use the A-type order, for which an altermagnetic
spin splitting is predicted. It is not the lowest collinear order in DFT: a spin-degenerate stripe order
(Mavani et al., Phys. Rev. Materials 10, 064408 (2026)) lies 14.3 meV per formula unit below it at these
settings (`dft/stripe/`). With the defects relaxed in both orders, the formation energies of V<sub>Mn</sub> and
Se<sub>Mn</sub> shift by +0.26 and -0.01 eV (`data/order_shift.csv`). The stripe-order V<sub>Mn</sub> converges to
M = -3.06 μ<sub>B</sub>; at its relaxed geometry M = -5 and -7 lie 0.21 and 0.45 eV higher (`dft/stripe/sm_*`). Magnetometry shows w-MnSe is antiferromagnetic below 62 K (Grzybowski et al.,
arXiv:2611.09509), but the type of order has not been determined.

**Authors:** Tomás Rojas Solorzano (CICIMA and Sede de Occidente, Universidad de Costa Rica;
tomas.rojas_s@ucr.ac.cr) and Dale Igram (University of Michigan–Dearborn).

**Paper:** T. Rojas Solorzano and D. Igram, *Native point defects in rock-salt and wurtzite MnSe: formation
energies across the growth window and spin splitting in the A-type order*, arXiv:2611.08986
(https://arxiv.org/abs/2611.08986). This repository is the data deposit for the version submitted to
Physica B: Condensed Matter; tagged releases correspond to manuscript versions.

> **License:** the data, calculation inputs and outputs, and figures are released under
> [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) (`LICENSE`); the scripts (`*.py`, `*.sh`) under the
> MIT License (`LICENSE-CODE`). If you use them, please cite the paper and this repository.

## What is in here

| Path | Contents |
|---|---|
| `data/total_energies.csv` | The raw total energies (Ry, 120 Ry cutoff) behind every formation energy, with the output file each one comes from |
| `data/formation_energies.csv` | Formation energies (eV) of the four defects in both polymorphs at four chemical-potential limits |
| `data/local_structure.csv` | First-shell relaxation, largest displacement, and Mn site-moment changes per defect |
| `data/order_shift.csv` | Change of the formation energies between the A-type and stripe orders, unrelaxed and relaxed |
| `data/relaxation_energies.csv` | Relaxation energy of every defect (ideal sites to relaxed, 70 Ry), with the two total energies and the files they come from, and the polymorph difference with and without relaxation. The wurtzite `gain_eV` excludes the 0.032 eV the 72-atom host gains from u = 0.375 to 0.378 (`dft/analysis/relax_gains.py`) |
| `data/structures/` | Relaxed defect supercells (extended XYZ, Å) with per-site magnetic moments |
| `data/tab_*.tex` | The same tables as LaTeX, including the unfolded spin splitting (`tab_stage5.tex`) |
| `dft/` | Quantum ESPRESSO inputs (and selected outputs) for every stage, plus the scripts that build the cells and analyse the results |
| `figures/` | Formation-energy diagram, bulk and defect-cell spin splitting of A-type w-MnSe, magnetic orders, structure renderings, CIF files |

Inside `dft/`:

| Folder | Stage |
|---|---|
| `stage1/` | Cutoff convergence and magnetic-ordering tests |
| `usweep/` | Bulk relaxations for U = 0–6 eV (sensitivity) |
| `hp/` | Linear-response (DFPT, `hp.x`) determination of U |
| `wzord/`, `sym/` | Wurtzite magnetic orders; symmetry checks of the spin splitting (α-MnSe null control, k-circle) |
| `bands/` | Bulk relaxations and the w-MnSe band structure along L₁–Γ–L₂ |
| `stage3prep/` | Elemental references (trigonal Se, Mn) and supercell preparation |
| `defects/`, `prod/`, `anti/` | Defect supercells: initial inputs, production relaxations (70 Ry) and single points (120 Ry), inputs and outputs, and the constrained-moment scans for Mn<sub>Se</sub> (`anti/`, `prod/wzscan/`) |
| `stage5/` | Spin splitting in supercells: band unfolding scripts, spectral weights and bands of the pristine and defect cells, peak tracking, the comparison of estimators (`estimator_compare.py`). The wavefunctions were deleted after analysis; `RECREATE_WAVEFUNCTIONS.md` gives the exact steps, cost and checks to regenerate them |
| `stripe/` | A-type, stripe and ferromagnetic orders of w-MnSe in a doubled cell; defects in both orders in a 64-atom cell, unrelaxed (`sd_*`), relaxed in legs (`sr_*_L*`, `run_relax.sh`, `order_shift_relaxed.py`) and spin-state checks (`sm_*`) |
| `spinscan/` | Constrained-magnetization scans at the relaxed defect geometries, and the off-center relaxation of rock-salt Se<sub>Mn</sub> |
| `unrelaxed/` | Single points of the 72-atom wurtzite V<sub>Se</sub> and Se<sub>Mn</sub> cells at ideal sites, the start of their relaxations |
| `prod_relaxedu/` | 72-atom wurtzite host at the relaxed internal parameter, the reference for the formation energies |
| `analysis/` | Per-site positions and moments of every defect cell, the local-structure analysis, and the relaxation energies (`relax_gains.py`, which writes `data/relaxation_energies.csv`) |

## Calculation settings (summary)

- Quantum ESPRESSO 7.5, PBE, PAW datasets from pslibrary (`Mn.pbe-spn-kjpaw_psl.0.3.1`,
  `Se.pbe-dn-kjpaw_psl.1.0.0`), spin-polarized, collinear, no spin–orbit coupling.
- Cutoffs 120 Ry / 960 Ry for all reported energies; defect relaxations at 70 Ry / 560 Ry.
- DFT+U (Dudarev), U = 3.7574 eV on Mn 3d from linear response, ortho-atomic projectors.
- α-MnSe in AFM-II order (64-atom 2×2×2 supercell of the conventional cell); w-MnSe in the
  A-type order (72-atom 3×3×2 supercell). Γ-centred 3×3×3 k-grids for the supercells.
- Chemical potentials: Se-rich from trigonal Se (also bounded by MnSe₂ formation); Mn-rich from
  the experimental formation enthalpy of MnSe, ΔH<sub>f</sub> = −1.70 eV.

## Reproducing the tables

```bash
cd dft
python3 formation_energies.py        # all 16 formation energies from the raw totals
python3 make_formation_table.py      # LaTeX table + lowest-energy defect vs chemical potential
python3 export_data.py               # regenerates everything in data/
python3 analysis/relax_gains.py ../data/relaxation_energies.csv   # relaxation energies, ideal-site differences
cd analysis && python3 local_analysis.py
```

Python 3 with NumPy (and Matplotlib for the figures). Paths to cluster directories inside the
Slurm scripts refer to the machine where the calculations ran.
