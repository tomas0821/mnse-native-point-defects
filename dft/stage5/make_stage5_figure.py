#!/usr/bin/env python3
"""Fig. — unfolded spin splitting of w-MnSe primitive band 27 along L1-Gamma-L2, by peak tracking:
the pristine curve at the relaxed internal parameter (s_relaxed_u.npz, written by
make_stage5_table.py, the like-for-like reference) vs the four defect cells (s_<defect>.npz from
peak_track.py --dump). Panel (b) spread: Gaussian widths 0.02-0.05 eV (spread_<defect>.npz). Writes
../../figures/fig_stage5_splitting.{pdf,png}."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "figures", "fig_stage5_splitting")
INK, INK2, GRID = "#1f1f1e", "#5f5e5a", "#e4e3df"
# Same colour + dash per defect as the formation-energy figure (validated palette slots 1-4);
# pristine is the neutral reference in ink grey.
STYLE = {"V_Mn": ("#2a78d6", "-"), "V_Se": ("#eb6834", (0, (5, 2))),
         "Se_Mn": ("#1baf7a", (0, (5, 1.5, 1, 1.5))), "Mn_Se": ("#eda100", (0, (1.2, 1.4)))}
LABEL = {"pristine": "pristine ($u = 0.378$)", "V_Mn": r"$V_{\mathrm{Mn}}$", "V_Se": r"$V_{\mathrm{Se}}$",
         "Se_Mn": r"$\mathrm{Se}_{\mathrm{Mn}}$", "Mn_Se": r"$\mathrm{Mn}_{\mathrm{Se}}$"}
plt.rcParams.update({"font.size": 8.5, "axes.edgecolor": INK2, "axes.labelcolor": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "font.family": "DejaVu Sans",
                     "mathtext.fontset": "dejavusans"})

x = np.concatenate([np.linspace(0, 1, 21), 1 + np.linspace(0, 1, 21)[1:]])  # L1..G..L2, 41 pts
fig, (ax, bx) = plt.subplots(1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1.6, 1], "wspace": 0.35})
ax.axhline(0, color=GRID, lw=0.8, zorder=0)
ax.axvspan(1 - 3.5 / 20, 1 + 3.5 / 20, color="#efeeea", lw=0, zorder=0)  # near-degenerate: not interpreted
ax.text(1.0, -285, "bands 27, 28\nclose", ha="center", va="bottom", fontsize=6.5, color=INK2)
present = []
for key in ("pristine", "V_Mn", "V_Se", "Se_Mn", "Mn_Se"):
    f = os.path.join(HERE, "s_relaxed_u.npz" if key == "pristine" else f"s_{key}.npz")
    if not os.path.exists(f):
        continue
    d = np.load(f)
    s = d["S"]
    if key == "pristine":
        ax.plot(x, s, color=INK2, lw=2.6, alpha=0.35, zorder=2, label=LABEL[key], solid_capstyle="round")
    else:
        c, ls = STYLE[key]
        ax.plot(x, s, color=c, ls=ls, lw=1.8, zorder=3, label=LABEL[key])
    present.append((key, float(np.abs(d["odd"][13])), float(np.abs(d["even"][13]))))  # at the pristine peak
ax.set_xticks([0, 1, 2], [r"L$_1$", r"$\Gamma$", r"L$_2$"])
ax.set_xlim(0, 2)
ax.set_ylabel(r"$E_\uparrow - E_\downarrow$ (meV)")
ax.set_ylim(-300, 300)
ax.legend(frameon=False, fontsize=7.5, loc="upper left", ncol=1)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
ax.text(0.0, 1.03, "(a)", transform=ax.transAxes, fontsize=8.5, color=INK)

# (b) altermagnetic (odd) vs uniform (even) parts, grouped bars per cell
keys = [k for k, _, _ in present]
XT = {"pristine": "pristine", **{k: LABEL[k] for k in STYLE}}
odd = [o for _, o, _ in present]
even = [e for _, _, e in present]
xi = np.arange(len(keys))
cols = [INK2 if k == "pristine" else STYLE[k][0] for k in keys]
bx.bar(xi - 0.19, odd, width=0.36, color=cols, edgecolor="white", linewidth=1.0, zorder=3)
bx.bar(xi + 0.19, even, width=0.36, color=cols, edgecolor="white", linewidth=1.0, hatch="////",
       alpha=0.55, zorder=3)
for i, k in enumerate(keys):
    g = os.path.join(HERE, f"spread_{k}.npz")
    if k == "pristine" or not os.path.exists(g):
        continue
    d = np.load(g)
    for off, a, val in ((-0.19, d["odd"], odd[i]), (0.19, d["even"], even[i])):
        bx.plot([i + off, i + off], [a.min(), a.max()], color=INK, lw=0.9, zorder=4)
        for y in (a.min(), a.max()):
            bx.plot([i + off - 0.06, i + off + 0.06], [y, y], color=INK, lw=0.9, zorder=4)
bx.set_xticks(xi, [XT[k] for k in keys], fontsize=7.5)
bx.set_ylabel("Splitting at the pristine peak (meV)")
bx.set_ylim(0, 340)
bx.grid(axis="y", color=GRID, lw=0.6, zorder=0)
for sp in ("top", "right"):
    bx.spines[sp].set_visible(False)
bx.text(0.0, 1.03, "(b)", transform=bx.transAxes, fontsize=8.5, color=INK)
from matplotlib.patches import Patch  # noqa: E402
bx.legend(handles=[Patch(facecolor=INK2, edgecolor="white", label="odd"),
                   Patch(facecolor=INK2, edgecolor="white", hatch="////", alpha=0.55, label="even")],
          frameon=False, fontsize=7, loc="upper right")
fig.subplots_adjust(left=0.095, right=0.98, bottom=0.14, top=0.9)
fig.savefig(OUT + ".pdf")
fig.savefig(OUT + ".png", dpi=300)
print("wrote", os.path.abspath(OUT), "| cells:", ", ".join(f"{k} odd {o:.1f} even {e:.1f}" for k, o, e in present))
