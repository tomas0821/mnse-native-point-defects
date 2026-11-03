#!/usr/bin/env python3
"""Neutral native-defect formation energies in alpha-MnSe (rock-salt, 64-atom) and w-MnSe
(wurtzite, 72-atom) from the raw 120 Ry totals. Single source for the manuscript table/figure.

E_f(dmu) = E_def - E_host - sum_i n_i mu_i ,  n_i = atoms ADDED (removed < 0)
mu_Se = mu_Se^0 + dmu (mu_Se^0 = trigonal Se / atom), mu_Mn = E_MnSe(host) - mu_Se
  =>  E_f(dmu) = E_f(0) + (n_Mn - n_Se) * dmu
Limits: Se-rich dmu = 0 (elemental Se; Aapro comparison); MnSe2-bounded Se-rich dmu = -0.158 eV
(expt [78Luk]) or -0.050 eV (our DFT); Mn-rich dmu = DH_f = -1.70 eV (expt 298 K, FERE-style
anchor, notebook 2026-08-12). Each polymorph is in equilibrium with ITSELF (own E_MnSe); the
same window is used for both (DFT polymorph ordering is unreliable, see notebook).
Energies (Ry) harvested 2026-11-03 from /scratch/global/tomas.rojas_s/mnse/{prod,anti,stage3prep}.
"""
RY = 13.605693123

E_SE3 = -1205.44004815                 # stage3prep/s3_ref_Se.out, trigonal Se, 3 atoms
HOST = {"rs": (-22026.89992696, 32),   # prod/d_rs_pristine_s120.out, 32 f.u.
        "wz": (-24780.54010775, 36)}   # prod/d_wz_pristine_s120.out, 36 f.u.
DEF = {  # (E_def Ry, n_Mn, n_Se, source, constraint/total magnetization)
    "rs": {"V_Mn": (-21740.31691645, -1, 0, "prod/d_rs_V_Mn_s120.out", "M=3.71"),
           "V_Se": (-21624.82558186, 0, -1, "prod/d_rs_V_Se_s120.out", "M=0"),
           "Se_Mn": (-22142.04600122, -1, 1, "prod/d_rs_Se_Mn_s120.out", "M=3.00"),
           "Mn_Se": (-21911.10496013, 1, -1, "anti/rsMnSe_M5_s120.out", "M=5 (constrained)")},
    "wz": {"V_Mn": (-24493.87426264, -1, 0, "prod/d_wz_V_Mn_s120.out", "M=-5.25"),
           "V_Se": (-24378.46008823, 0, -1, "prod/d_wz_V_Se_s120.out", "M=0"),
           "Se_Mn": (-24895.70752742, -1, 1, "prod/d_wz_Se_Mn_s120.out", "M=-5.00"),
           "Mn_Se": (-24664.70559843, 1, -1, "prod/d_wz_Mn_Se_s120.out", "M=3 (constrained)")},
}
LIMITS = {"Se-rich (elem. Se)": 0.0, "Se-rich (MnSe2, expt)": -0.158,
          "Se-rich (MnSe2, DFT)": -0.050, "Mn-rich": -1.70}
DMU_MNRICH = -1.70


def table():
    mu_se = E_SE3 / 3
    out = {}
    for p, (eh, nfu) in HOST.items():
        mu_mn = eh / nfu - mu_se
        for d, (ed, nmn, nse, _, _) in DEF[p].items():
            ef0 = (ed - eh - nmn * mu_mn - nse * mu_se) * RY
            out[(p, d)] = (ef0, nmn - nse)
    return out


def ef(t, p, d, dmu):
    ef0, slope = t[(p, d)]
    return ef0 + slope * dmu


if __name__ == "__main__":
    t = table()
    print(f"mu_Se = {E_SE3/3:.8f} Ry")
    hdr = "".join(f"{k:>24s}" for k in LIMITS)
    print(f"{'defect':10s}{'slope':>6s}{hdr}")
    for p in ("rs", "wz"):
        for d in ("V_Mn", "V_Se", "Se_Mn", "Mn_Se"):
            row = "".join(f"{ef(t, p, d, v):24.3f}" for v in LIMITS.values())
            print(f"{p}:{d:7s}{t[(p, d)][1]:6d}{row}")
    print("wz - rs, Se-rich: " + ", ".join(
        f"{d} {ef(t,'wz',d,0)-ef(t,'rs',d,0):+.3f}" for d in ("V_Mn", "V_Se", "Se_Mn", "Mn_Se")))
