#!/usr/bin/env python3
"""
Structure figures for the lab notebook. Coordinates come from the SAME builders
used to generate the QE inputs (dft/gen_stage1.py), so the figures show what was
actually computed, not an idealised sketch.

Palette validated for colour-vision deficiency (see notebook): normal-vision
pairwise dE_OKLab >= 15 and worst-case dichromat dE >= 8. Spin is ALSO encoded by
an arrow, so identity never rests on colour alone.

    python3 dft/make_structure_figures.py
"""
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(os.path.dirname(HERE), "figures")
os.makedirs(FIGS, exist_ok=True)

MN_UP, MN_DN, SE = "#2E5FA3", "#C1442E", "#A8AEB8"
INK, MUTED = "#1F2328", "#6B7280"

A_RS = 5.4616
A_WZ, C_WZ, U_WZ = 4.19, 6.68, 0.375


# ---------- structures (mirrors gen_stage1.py) ----------
def rocksalt(n, order):
    L, h = 2 * n, A_RS / 2
    mn, se = [], []
    for i in range(L + 1):
        for j in range(L + 1):
            for k in range(L + 1):
                if (i + j + k) % 2 == 0:
                    s = (+1 if order == "FM" else
                         (-1 if k % 2 else +1) if order == "AFMI" else
                         (-1 if ((i + j + k) // 2) % 2 else +1))
                    mn.append((np.array([i, j, k]) * h, s))
                else:
                    se.append(np.array([i, j, k]) * h)
    return mn, se, A_RS * n


def wurtzite(nx=2, ny=2, nz=1):
    a1 = np.array([A_WZ, 0, 0])
    a2 = np.array([-A_WZ / 2, A_WZ * math.sqrt(3) / 2, 0])
    a3 = np.array([0, 0, C_WZ])
    mn_f = [(np.array([1 / 3, 2 / 3, 0.0]), +1), (np.array([2 / 3, 1 / 3, 0.5]), -1)]
    se_f = [np.array([1 / 3, 2 / 3, U_WZ]), np.array([2 / 3, 1 / 3, U_WZ + 0.5])]
    mn, se = [], []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):   # range(nz), not nz+1 -- else a whole layer sits
                                  # outside the drawn cell and looks detached
                t = i * a1 + j * a2 + k * a3
                for f, s in mn_f:
                    mn.append((f[0] * a1 + f[1] * a2 + f[2] * a3 + t, s))
                for f in se_f:
                    se.append(f[0] * a1 + f[1] * a2 + f[2] * a3 + t)
    return mn, se, (a1, a2, a3)


# ---------- rendering ----------
def view_dir(elev, azim):
    e, a = math.radians(elev), math.radians(azim)
    return np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])


def draw(ax, mn, se, elev, azim, r_mn=170, r_se=90, arrows=True, se_alpha=0.85):
    """Painter's algorithm: plot far atoms first so near ones overlap correctly."""
    v = view_dir(elev, azim)
    items = [(p, "se", None) for p in se] + [(p, "mn", s) for p, s in mn]
    items.sort(key=lambda it: float(np.dot(it[0], v)))
    for p, kind, s in items:
        if kind == "se":
            ax.scatter(*p, s=r_se, c=SE, edgecolors="white", linewidths=0.6,
                       alpha=se_alpha, depthshade=False, zorder=2)
        else:
            ax.scatter(*p, s=r_mn, c=(MN_UP if s > 0 else MN_DN),
                       edgecolors="white", linewidths=0.8, depthshade=False, zorder=3)
    if arrows:
        for p, s in mn:
            ax.quiver(p[0], p[1], p[2] - 0.55 * s, 0, 0, 1.1 * s,
                      color=INK, linewidth=1.0, arrow_length_ratio=0.45, zorder=4)


def cube_edges(ax, L):
    pts = np.array([[0, 0, 0], [L, 0, 0], [L, L, 0], [0, L, 0],
                    [0, 0, L], [L, 0, L], [L, L, L], [0, L, L]], float)
    E = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
         (0, 4), (1, 5), (2, 6), (3, 7)]
    for i, j in E:
        ax.plot(*zip(pts[i], pts[j]), color=MUTED, lw=0.8, alpha=0.55, zorder=1)


def hex_edges(ax, vecs, nx, ny, nz):
    """Outline the FULL nx x ny x nz block, not just one primitive cell."""
    a1, a2, a3 = vecs
    base = [np.zeros(3), nx * a1, nx * a1 + ny * a2, ny * a2]
    for k in (0, nz):
        loop = [b + k * a3 for b in base] + [base[0] + k * a3]
        for p, q in zip(loop, loop[1:]):
            ax.plot(*zip(p, q), color=MUTED, lw=0.8, alpha=0.55, zorder=1)
    for b in base:
        ax.plot(*zip(b, b + nz * a3), color=MUTED, lw=0.8, alpha=0.55, zorder=1)


def tidy(ax, pts, elev, azim, title, sub, aspect=None):
    P = np.array(pts)
    lo, hi = P.min(0), P.max(0)
    c, ext = (lo + hi) / 2, (hi - lo)
    pad = ext.max() * 0.08
    ax.set_xlim(c[0] - ext[0] / 2 - pad, c[0] + ext[0] / 2 + pad)
    ax.set_ylim(c[1] - ext[1] / 2 - pad, c[1] + ext[1] / 2 + pad)
    ax.set_zlim(c[2] - ext[2] / 2 - pad, c[2] + ext[2] / 2 + pad)
    # box aspect must follow the real extents, or a tall cell wastes the canvas
    ax.set_box_aspect(aspect if aspect else tuple(np.maximum(ext, ext.max() * 0.35)))
    ax.view_init(elev=elev, azim=azim)
    ax.set_axis_off()
    ax.set_title(title, fontsize=12, color=INK, pad=2, weight="bold")
    if sub:
        ax.text2D(0.5, 0.955, sub, transform=ax.transAxes, ha="center",
                  fontsize=9, color=MUTED)


LEGEND = [Line2D([], [], marker="o", ls="", ms=9, mfc=MN_UP, mec="white", label="Mn  spin up (↑)"),
          Line2D([], [], marker="o", ls="", ms=9, mfc=MN_DN, mec="white", label="Mn  spin down (↓)"),
          Line2D([], [], marker="o", ls="", ms=7, mfc=SE, mec="white", label="Se")]


def fig_rocksalt():
    mn, se, L = rocksalt(1, "AFMII")
    fig = plt.figure(figsize=(6.4, 6.0), dpi=300)
    ax = fig.add_subplot(111, projection="3d")
    draw(ax, mn, se, 18, -60)
    cube_edges(ax, L)
    tidy(ax, [p for p, _ in mn] + se, 18, -60,
         "α-MnSe — rock-salt (Fm$\\bar{3}$m), AFM-II",
         f"a = {A_RS} Å (Pradhan 2025, XRD)")
    ax.text2D(0.5, 0.045,
              "AFM-II spin pattern drawn on the conventional cubic cell for clarity.\n"
              "That cell is NOT a valid magnetic cell — type-II order needs the [111]\n"
              "doubling; the calculations use the 2×2×2 (64-atom) supercell.",
              transform=ax.transAxes, ha="center", fontsize=7.8, color=MUTED)
    ax.legend(handles=LEGEND, loc="lower center", ncol=3, frameon=False,
              fontsize=8.5, bbox_to_anchor=(0.5, -0.05))
    fig.tight_layout()
    p = os.path.join(FIGS, "fig1_alpha-MnSe_rocksalt.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def fig_orderings():
    """Mn sublattice only -- Se removed so the magnetic pattern is unmistakable."""
    fig = plt.figure(figsize=(11.5, 5.8), dpi=300)
    for idx, (order, name, plane, note) in enumerate([
            ("AFMI", "AFM-I", "(001)", "alternating (001) planes  ·  4 same / 8 opposite NN"),
            ("AFMII", "AFM-II", "(111)", "alternating (111) planes  ·  6 same / 6 opposite NN")]):
        mn, se, L = rocksalt(1, order)
        ax = fig.add_subplot(1, 2, idx + 1, projection="3d")
        draw(ax, mn, [], 20, -55, r_mn=210)
        cube_edges(ax, L)
        gs = "  ← ground state (this work, +0 meV)" if order == "AFMII" else "  (+15.9 meV/f.u.)"
        tidy(ax, [p for p, _ in mn], 20, -55,
             f"{name} — {plane} stacking{gs}", note)
    fig.suptitle("α-MnSe magnetic orderings — Mn sublattice only (Se omitted for clarity)",
                 fontsize=12.5, color=INK, weight="bold", y=0.99)
    fig.legend(handles=LEGEND[:2], loc="lower center", ncol=2, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0.045, 1, 0.96))
    p = os.path.join(FIGS, "fig2_alpha-MnSe_orderings.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


def fig_wurtzite():
    mn, se, vecs = wurtzite(2, 2, 2)
    fig = plt.figure(figsize=(6.0, 7.0), dpi=300)
    ax = fig.add_subplot(111, projection="3d")
    draw(ax, mn, se, 14, -68)
    hex_edges(ax, vecs, 2, 2, 2)
    tidy(ax, [p for p, _ in mn] + se, 14, -68,
         "w-MnSe — wurtzite (P6$_3$mc), altermagnetic",
         f"a = {A_WZ} Å, c = {C_WZ} Å (Bezzerga 2025)  ·  2×2×2 cells")
    ax.text2D(0.5, 0.035,
              "the two Mn sublattices are related by a 6-fold rotation about c\n"
              "combined with a non-primitive translation — the altermagnetic symmetry",
              transform=ax.transAxes, ha="center", fontsize=8.2, color=MUTED)
    ax.legend(handles=LEGEND, loc="lower center", ncol=3, frameon=False,
              fontsize=8.5, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout()
    p = os.path.join(FIGS, "fig3_w-MnSe_wurtzite.png")
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return p


if __name__ == "__main__":
    for f in (fig_rocksalt(), fig_orderings(), fig_wurtzite()):
        print("wrote", f, f"({os.path.getsize(f)//1024} kB)")
