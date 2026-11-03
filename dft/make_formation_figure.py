#!/usr/bin/env python3
"""Fig. — neutral native-defect formation energies across the chemical-potential window,
alpha-MnSe vs w-MnSe. Numbers come only from formation_energies.py (raw 120 Ry totals).
Writes figures/fig_formation_energies.{pdf,png}."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from formation_energies import table, ef, DMU_MNRICH  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "fig_formation_energies")

# Categorical slots 1-4 of the validated reference palette (light; validator: CVD worst
# adjacent dE 9.1, normal 22.9; aqua/yellow < 3:1 contrast -> every line is direct-labelled
# and carries a distinct line style as secondary encoding).
STYLE = {"V_Mn": ("#2a78d6", "-"), "V_Se": ("#eb6834", (0, (5, 2))),
         "Se_Mn": ("#1baf7a", (0, (5, 1.5, 1, 1.5))), "Mn_Se": ("#eda100", (0, (1.2, 1.4)))}
LABEL = {"V_Mn": r"$V_{\mathrm{Mn}}$", "V_Se": r"$V_{\mathrm{Se}}$",
         "Se_Mn": r"$\mathrm{Se}_{\mathrm{Mn}}$", "Mn_Se": r"$\mathrm{Mn}_{\mathrm{Se}}$"}
INK, INK2, GRID = "#1f1f1e", "#5f5e5a", "#e4e3df"
DMU_MNSE2 = -0.158

plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "font.family": "DejaVu Sans",
                     "mathtext.fontset": "dejavusans"})

t = table()
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True,
                         gridspec_kw={"wspace": 0.30})
titles = {"rs": r"$\alpha$-MnSe (rock-salt, AFM-II)", "wz": "w-MnSe (wurtzite, altermagnet)"}
x0, x1 = DMU_MNRICH, 0.0
for ax, p in zip(axes, ("rs", "wz")):
    ax.axvspan(DMU_MNSE2, x1, color="#efeeea", lw=0, zorder=0)
    ax.text((DMU_MNSE2 + x1) / 2, 5.0, "MnSe$_2$\nforms", ha="center", va="center",
            fontsize=6.5, color=INK2)
    ax.grid(axis="y", color=GRID, lw=0.6, zorder=0)
    ends = []
    for d in ("V_Mn", "V_Se", "Se_Mn", "Mn_Se"):
        c, ls = STYLE[d]
        y0, y1 = ef(t, p, d, x0), ef(t, p, d, x1)
        ax.plot([x0, x1], [y0, y1], color=c, ls=ls, lw=2.0, solid_capstyle="round",
                zorder=3, label=LABEL[d])
        ends.append((y1, d))
    # direct labels at the Se-rich end, nudged apart if closer than 0.32 eV
    ends.sort()
    placed = []
    for y, d in ends:
        if placed and y - placed[-1] < 0.40:
            y = placed[-1] + 0.40
        placed.append(y)
        ax.text(x1 + 0.03, y, LABEL[d], va="center", ha="left", fontsize=8, color=INK)
    ax.set_xlim(x0, x1)
    ax.set_ylim(0, 8)
    ax.set_title(titles[p], fontsize=8.5, color=INK, loc="left")
    ax.set_xlabel(r"$\Delta\mu_{\mathrm{Se}}$ (eV)")
    ax.text(0.01, -0.2, "Mn-rich", transform=ax.transAxes, fontsize=7, color=INK2, ha="left")
    ax.text(0.99, -0.2, "Se-rich", transform=ax.transAxes, fontsize=7, color=INK2, ha="right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
axes[0].set_ylabel("Formation energy (eV)")
axes[0].legend(loc="upper left", frameon=False, fontsize=7.5, ncol=2, handlelength=2.6,
               bbox_to_anchor=(0.0, 0.93))
fig.subplots_adjust(left=0.08, right=0.92, bottom=0.2, top=0.9)
os.makedirs(os.path.dirname(OUT), exist_ok=True)
fig.savefig(OUT + ".pdf")
fig.savefig(OUT + ".png", dpi=300)
print("wrote", os.path.abspath(OUT) + ".{pdf,png}")
