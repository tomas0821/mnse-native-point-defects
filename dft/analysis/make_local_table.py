#!/usr/bin/env python3
"""Manuscript Table 2 (local relaxation and net magnetization) from local_analysis.py.
Net cell magnetizations |M| (mu_B) are the 'total magnetization' lines of the 120 Ry outputs
(prod/d_*_s120.out, anti/rsMnSe_M5_s120.out), harvested 2026-11-03; Mn_Se values are the
constrained minima of the M-scans. Writes ../../data/tab_local.tex."""
import os
import local_analysis as L

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "data", "tab_local.tex")
M = {("rs", "V_Mn"): "3.71", ("rs", "V_Se"): "0.00", ("rs", "Se_Mn"): "3.00", ("rs", "Mn_Se"): "5$^a$",
     ("wz", "V_Mn"): "5.25", ("wz", "V_Se"): "0.00", ("wz", "Se_Mn"): "5.00", ("wz", "Mn_Se"): "3$^a$"}
TEX = {"V_Mn": r"$V_{\mathrm{Mn}}$", "V_Se": r"$V_{\mathrm{Se}}$",
       "Se_Mn": r"$\mathrm{Se}_{\mathrm{Mn}}$", "Mn_Se": r"$\mathrm{Mn}_{\mathrm{Se}}$"}
rows = []
for d in ("V_Mn", "V_Se", "Se_Mn", "Mn_Se"):
    a, w = L.analyse("rs", d), L.analyse("wz", d)
    rows.append(f"{TEX[d]} & ${a['radial']:+.3f}$ & {a['maxdisp']:.3f} & {M[('rs', d)]} & "
                f"${w['radial']:+.3f}$ & {w['maxdisp']:.3f} & {M[('wz', d)]} \\\\")
tex = r"""\begin{table}[t]
\centering
\caption{Local relaxation and net magnetization of the defect cells. $\Delta r_1$: mean radial
displacement of the first coordination shell (positive outward); $|u|_{\max}$: largest atomic
displacement in the cell; $|M|$: net magnetization of the cell ($\mu_{\mathrm{B}}$).
First shells: six neighbors at 2.77~\AA\ in $\alpha$-MnSe, four at 2.57--2.61~\AA\ in w-MnSe.
$^a$Constrained at the minimum of the $M$ scan.}
\label{tab:local}
\begin{tabular}{lrrrrrr}
\toprule
 & \multicolumn{3}{c}{$\alpha$-MnSe} & \multicolumn{3}{c}{w-MnSe} \\
\cmidrule(lr){2-4}\cmidrule(lr){5-7}
Defect & $\Delta r_1$ (\AA) & $|u|_{\max}$ (\AA) & $|M|$ & $\Delta r_1$ (\AA) & $|u|_{\max}$ (\AA) & $|M|$ \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
open(OUT, "w").write(tex)
print(tex[tex.index("midrule"):tex.index("bottomrule")])
