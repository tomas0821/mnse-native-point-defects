#!/usr/bin/env python3
"""Manuscript Table 1 + dominant-defect crossovers, from formation_energies.py only.
Writes data/tab_formation.tex and prints the lowest-energy defect regimes."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from formation_energies import table, ef, DMU_MNRICH  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "tab_formation.tex")
D = ("V_Mn", "V_Se", "Se_Mn", "Mn_Se")
TEX = {"V_Mn": r"$V_{\mathrm{Mn}}$", "V_Se": r"$V_{\mathrm{Se}}$",
       "Se_Mn": r"$\mathrm{Se}_{\mathrm{Mn}}$", "Mn_Se": r"$\mathrm{Mn}_{\mathrm{Se}}$"}
LIM = [("elem", 0.0), ("MnSe2", -0.158), ("Mnrich", DMU_MNRICH)]
t = table()

rows = []
for d in D:
    a = [ef(t, "rs", d, v) for _, v in LIM]
    w = [ef(t, "wz", d, v) for _, v in LIM]
    dn = t[("rs", d)][1]
    rows.append(f"{TEX[d]} & ${dn:+d}$ & " + " & ".join(f"{x:.2f}" for x in a) + " & "
                + " & ".join(f"{x:.2f}" for x in w) + f" & ${w[0]-a[0]:+.2f}$ \\\\")
tex = r"""\begin{table*}[t]
\centering
\caption{Formation energies (eV) of the neutral native defects in $\alpha$-MnSe (64-atom cell)
and w-MnSe (72-atom cell). Se-rich: elemental trigonal Se ($\Delta\mu_{\mathrm{Se}} = 0$) and the
MnSe$_2$-bounded limit ($\Delta\mu_{\mathrm{Se}} = -0.158$~eV); Mn-rich:
$\Delta\mu_{\mathrm{Se}} = \Delta H_f = -1.70$~eV. The slope
$n_{\mathrm{Mn}} - n_{\mathrm{Se}}$ of $E_f$ versus $\Delta\mu_{\mathrm{Se}}$ is listed, and the
last column gives the polymorph difference $E_f(\mathrm{w}) - E_f(\alpha)$, which is the same in
every limit.}
\label{tab:Ef}
\begin{tabular}{lrrrrrrrr}
\toprule
 & & \multicolumn{3}{c}{$\alpha$-MnSe} & \multicolumn{3}{c}{w-MnSe} & \\
\cmidrule(lr){3-5}\cmidrule(lr){6-8}
Defect & $n_{\mathrm{Mn}}-n_{\mathrm{Se}}$ & Se & MnSe$_2$ & Mn-rich & Se & MnSe$_2$ & Mn-rich & w$-\alpha$ \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table*}
"""
open(OUT, "w").write(tex)
print("wrote", os.path.abspath(OUT))

# lowest-energy defect along the window, fine scan
import numpy as np
for p in ("rs", "wz"):
    grid = np.linspace(0.0, DMU_MNRICH, 17001)
    low = [min(D, key=lambda d: ef(t, p, d, x)) for x in grid]
    print(f"\n{p}: lowest-energy defect regimes (Se-rich -> Mn-rich)")
    start = 0
    for i in range(1, len(grid) + 1):
        if i == len(grid) or low[i] != low[start]:
            print(f"   {low[start]:6s}  dmu_Se from {grid[start]:+.3f} to {grid[i-1]:+.3f} eV"
                  f"   (E_f {ef(t,p,low[start],grid[start]):.2f} -> {ef(t,p,low[start],grid[i-1]):.2f})")
            start = i
