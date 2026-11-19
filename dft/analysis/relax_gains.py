#!/usr/bin/env python3
"""Relaxation energies of the four neutral defects in both polymorphs (70/560 Ry, production cells), and the
polymorph difference E_f(w) - E_f(alpha) with and without relaxation.

gain = E(ideal sites) - E(relaxed). The wurtzite defect cells start from ideal u = 0.375 while the host reference is
relaxed (u = 0.3778); the host's own u relaxation, 0.9 meV/f.u. x 36 f.u. = 0.032 eV, is subtracted from every
wurtzite gain so that only the defect's relaxation is counted.
Unrelaxed difference = relaxed difference (Table 1, 120 Ry) + gain_w - gain_alpha.
usage: relax_gains.py [out.csv]   (the CSV lists the two energies behind every gain and the files they come from)"""
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from formation_energies import table, ef  # noqa: E402

RY = 13.605693123
U_HOST_W = 0.032   # eV, ideal-u penalty of the 72-atom wurtzite host


def energies(*files):
    E = []
    for f in files:
        E += [float(x) for x in re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)",
                                            open(os.path.join(HERE, "..", f), errors="ignore").read(), re.M)]
    return E


def gain(ideal, relaxed):
    """ideal, relaxed: tuples of output files; first energy of `ideal`, last energy of `relaxed`."""
    return (energies(*ideal)[0] - energies(*relaxed)[-1]) * RY


SRC = {  # (ideal-site run(s), relaxed run(s))
    ("rs", "V_Mn"): (("prod/d_rs_V_Mn_r70.out",), ("prod/d_rs_V_Mn_r70.out",)),
    ("rs", "V_Se"): (("prod/d_rs_V_Se_r70.out",), ("prod/d_rs_V_Se_r70.out",)),
    ("rs", "Se_Mn"): (("prod/d_rs_Se_Mn_r70.out",), ("prod/d_rs_Se_Mn_r70.out",)),
    ("rs", "Mn_Se"): (("anti/anti_M5p0.out",), ("anti/rsMnSe_M5.out",)),
    ("wz", "V_Mn"): (("prod/d_wz_V_Mn_r70.out",), ("prod/d_wz_V_Mn_r70.out", "prod/d_wz_V_Mn_rs.out")),
    ("wz", "V_Se"): (("unrelaxed/u_wz_V_Se.out",), ("prod/d_wz_V_Se_r70.out", "prod/d_wz_V_Se_rs.out")),
    ("wz", "Se_Mn"): (("unrelaxed/u_wz_Se_Mn.out",), ("prod/d_wz_Se_Mn_r70.out", "prod/d_wz_Se_Mn_rs.out")),
    ("wz", "Mn_Se"): (("prod/wzscan/wzanti_M3.out",), ("prod/wzMnSe_M3rs.out",)),
}

if __name__ == "__main__":
    import csv
    t = table()
    print("defect  gain_alpha  gain_w   w-alpha relaxed  w-alpha unrelaxed   (eV)")
    out = []
    for d in ("V_Mn", "V_Se", "Se_Mn", "Mn_Se"):
        ga = gain(*SRC[("rs", d)])
        gw = gain(*SRC[("wz", d)]) - U_HOST_W
        rel = ef(t, "wz", d, 0.0) - ef(t, "rs", d, 0.0)
        print(f"{d:6s}  {ga:10.3f}  {gw:7.3f}  {rel:+15.3f}  {rel + gw - ga:+17.3f}")
        for p, g in (("rs", ga), ("wz", gw)):
            ideal, relaxed = SRC[(p, d)]
            out.append([p, d, f"{energies(*ideal)[0]:.8f}", f"{energies(*relaxed)[-1]:.8f}", f"{g:.3f}",
                        f"{rel:+.3f}", f"{rel + gw - ga:+.3f}", " ".join(ideal), " ".join(relaxed)])
    if len(sys.argv) > 1:
        with open(sys.argv[1], "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["polymorph", "defect", "E_ideal_Ry", "E_relaxed_Ry", "gain_eV",
                        "w_minus_alpha_relaxed_eV", "w_minus_alpha_ideal_sites_eV", "ideal_from", "relaxed_from"])
            w.writerows(out)
        print("wrote", sys.argv[1])
