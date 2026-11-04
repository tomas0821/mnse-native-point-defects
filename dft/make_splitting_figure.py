#!/usr/bin/env python3
"""
The headline validation figure: w-MnSe altermagnetic spin splitting.

(a) spin-resolved bands along L1-Gamma-L2
(b) splitting vs azimuth on a k-circle: the sign alternates every 60 degrees

Same spin colours as the structure figures (validated for colour-vision
deficiency: normal-vision dE_OKLab 27.4, worst-case dichromat 18.2). Spin is
also encoded by line style and by explicit labels, never colour alone.

    python3 dft/make_splitting_figure.py
"""
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(os.path.dirname(HERE), "figures")

UP, DN = "#2E5FA3", "#C1442E"
INK, MUTED, GRID = "#1F2328", "#6B7280", "#D9DDE3"
NV = 31          # 62 valence electrons -> 31 occupied bands per spin channel


def read(path):
    t = open(path, errors="ignore").read()
    up = t.split("SPIN UP")[1].split("SPIN DOWN")[0]
    dn = t.split("SPIN DOWN")[1]
    grab = lambda b: [[float(x) for x in m.group(1).split()]
                      for m in re.finditer(r"bands \(ev\):\s*\n\n(.*?)(?=\n\s*\n)", b, re.S)]
    return np.array(grab(up)), np.array(grab(dn))


U, D = read(os.path.join(HERE, "bands", "b_wz_bands.out"))
Ug, Dg = read(os.path.join(HERE, "sym", "gw_bands.out"))

fig = plt.figure(figsize=(11.2, 4.9), dpi=300)

# ---------------- (a) bands along L1 - Gamma - L2 ----------------
ax = fig.add_subplot(1, 2, 1)
n = len(U)
x = np.arange(n)
vbm = max(U[:, NV - 1].max(), D[:, NV - 1].max())
lo, hi = NV - 7, NV
for b in range(lo, hi):
    ax.plot(x, U[:, b] - vbm, color=UP, lw=1.4, zorder=3,
            label="spin ↑" if b == lo else None)
    ax.plot(x, D[:, b] - vbm, color=DN, lw=1.4, ls=(0, (4, 2)), zorder=3,
            label="spin ↓" if b == lo else None)
# shade where the two channels separate, on the band carrying the effect
band = 26
ax.fill_between(x, U[:, band] - vbm, D[:, band] - vbm, color="#8A8F98",
                alpha=0.45, lw=0, zorder=2)
_w = np.concatenate([U[:, lo:hi], D[:, lo:hi]]) - vbm
ax.set_ylim(_w.min() - 0.25, _w.max() + 0.35)
for xt in (0, n // 2, n - 1):
    ax.axvline(xt, color=GRID, lw=0.9, zorder=1)
ax.set_xticks([0, n // 2, n - 1])
ax.set_xticklabels(["L$_1$", "Γ", "L$_2$"], fontsize=11)
ax.set_xlim(0, n - 1)
ax.set_ylabel("E − E$_{VBM}$  (eV)", fontsize=10, color=INK)
ax.set_title("(a)", fontsize=11, color=INK, loc="left")
ax.tick_params(labelsize=9, colors=MUTED)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax.spines[s].set_color(GRID)
ax.legend(frameon=False, fontsize=9, loc="upper center", ncol=2,
          bbox_to_anchor=(0.5, 1.13))

# ---------------- (b) azimuthal scan: the 6-fold sign alternation ------------
ax2 = fig.add_subplot(1, 2, 2)
m = len(Ug)
phi = np.linspace(0, 360, m)
spl = (Ug[:, band] - Dg[:, band]) * 1000
ax2.axhline(0, color=GRID, lw=1.0, zorder=1)
ax2.fill_between(phi, 0, spl, where=spl >= 0, color=UP, alpha=0.30, lw=0, zorder=2)
ax2.fill_between(phi, 0, spl, where=spl <= 0, color=DN, alpha=0.30, lw=0, zorder=2)
ax2.plot(phi, spl, color=INK, lw=1.6, zorder=4)
for z in range(0, 361, 60):
    ax2.axvline(z, color=GRID, lw=0.9, ls=(0, (2, 3)), zorder=1)
ax2.set_xticks(range(0, 361, 60))
ax2.set_xlim(0, 360)
ax2.set_xlabel("azimuth φ around the c-axis  (degrees)", fontsize=10, color=INK)
ax2.set_ylabel("E$_↑$ − E$_↓$   (meV)", fontsize=10, color=INK)
ax2.set_title("(b)", fontsize=11, color=INK, loc="left")
ax2.tick_params(labelsize=9, colors=MUTED)
for s in ("top", "right"):
    ax2.spines[s].set_visible(False)
for s in ("left", "bottom"):
    ax2.spines[s].set_color(GRID)
pk = np.abs(spl).max()
fig.tight_layout()
p = os.path.join(FIGS, "fig4_w-MnSe_spin-splitting.png")
fig.savefig(p, bbox_inches="tight", facecolor="white")
print("wrote", p, f"({os.path.getsize(p)//1024} kB)")
print(f"  peak splitting on path   : {np.abs((U[:,band]-D[:,band])*1000).max():.1f} meV")
print(f"  peak splitting on circle : {pk:.1f} meV")
