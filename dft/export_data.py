#!/usr/bin/env python3
"""Write the summary data files in ../data/ from the single-source scripts:
  total_energies.csv        raw 120 Ry total energies used for every formation energy
  formation_energies.csv    16 neutral formation energies x 4 chemical-potential limits
  local_structure.csv       first-shell relaxation, largest displacement, site moments
  structures/*.xyz          relaxed defect supercells (extended XYZ, Angstrom, with moments)
usage: python3 export_data.py      (run from dft/)"""
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "analysis"))
import formation_energies as F  # noqa: E402

DATA = os.path.join(HERE, "..", "data")
os.makedirs(os.path.join(DATA, "structures"), exist_ok=True)
DEF = ("V_Mn", "V_Se", "Se_Mn", "Mn_Se")
POL = {"rs": "alpha-MnSe (rock salt, 64-atom)", "wz": "w-MnSe (wurtzite, 72-atom)"}

with open(os.path.join(DATA, "total_energies.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["system", "polymorph", "E_total_Ry", "n_formula_units", "source_output", "net_magnetization"])
    w.writerow(["trigonal Se (3 atoms)", "-", F.E_SE3, "", "stage3prep/s3_ref_Se.out", ""])
    for p, (eh, nfu) in F.HOST.items():
        w.writerow(["pristine host", p, eh, nfu, f"prod/d_{p}_pristine_s120.out", "0"])
        for d in DEF:
            ed, nmn, nse, src, mag = F.DEF[p][d]
            w.writerow([d, p, ed, "", src, mag])

t = F.table()
with open(os.path.join(DATA, "formation_energies.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["polymorph", "defect", "slope_nMn_minus_nSe"] + [f"Ef_eV[{k}; dmu_Se={v}]" for k, v in F.LIMITS.items()])
    for p in ("rs", "wz"):
        for d in DEF:
            w.writerow([p, d, t[(p, d)][1]] + [f"{F.ef(t, p, d, v):.3f}" for v in F.LIMITS.values()])

cwd = os.getcwd()
os.chdir(os.path.join(HERE, "analysis"))
import local_analysis as L  # noqa: E402
with open(os.path.join(DATA, "local_structure.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["polymorph", "defect", "first_shell_n", "first_shell_species", "first_shell_R_A",
                "first_shell_mean_radial_disp_A", "max_atomic_disp_A", "bulk_site_moment_muB",
                "max_abs_moment_change_excl_antisite_muB", "at_R_A", "antisite_Mn_moment_muB"])
    for p in ("rs", "wz"):
        for d in DEF:
            o = L.analyse(p, d)
            w.writerow([p, d, o["shell_n"], o["shell_sp"], f"{o['shell_R']:.3f}", f"{o['radial']:+.3f}",
                        f"{o['maxdisp']:.3f}", f"{o['ref']:.3f}", f"{o['max_dm']:.3f}", f"{o['max_dm_R']:.2f}",
                        f"{o['antisite_m']:+.3f}" if "antisite_m" in o else ""])
            D = L.load(f"d_{p}_{d}.dat")
            c = D["CELL"]
            lat = " ".join(f"{x:.8f}" for x in c.flatten())
            with open(os.path.join(DATA, "structures", f"{p}_{d}_relaxed.xyz"), "w") as g:
                g.write(f"{len(D['FINAL'])}\n")
                g.write(f'Lattice="{lat}" Properties=species:S:1:pos:R:3:magmom:R:1 '
                        f'pbc="T T T" system="{POL[p]} {d}, relaxed geometry; site moments as in analysis/d_{p}_{d}.dat"\n')
                for (sp, r), m in zip(D["FINAL"], D["MAGN"]):
                    g.write(f"{sp:2s} {r[0]:14.8f} {r[1]:14.8f} {r[2]:14.8f} {m:9.4f}\n")
os.chdir(cwd)
print("wrote", sorted(os.listdir(DATA)), "and", len(os.listdir(os.path.join(DATA, "structures"))), "structures")
