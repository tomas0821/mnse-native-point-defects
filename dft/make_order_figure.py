#!/usr/bin/env python3
"""Fig.: the wurtzite lattice and the two magnetic orders compared in the paper.
(a) side view with the Mn sublattices of the A-type order; (b) A-type and (c) stripe AFMs-alpha, seen along
c, Mn only (layer A large, layer B small). Geometry: relaxed cell of usweep/wz_U37574.out, as in gen_stripe.py.
Writes ../figures/fig_orders.{pdf,png}."""
import itertools
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "fig_orders")
UP, DN, SE, INK, MUTED = "#2E5FA3", "#C1442E", "#A9AEB6", "#1F2328", "#6B7280"
A = np.array([[4.270583028, 0, 0], [-2.135291514, 3.698433391, 0], [0, 0, 6.854295104]])
MN = np.array([[0.0, 2.4656222608, -0.0096777392], [2.1352915140, 1.2328111304, 3.4174697914]])
SEL = np.array([[0.0, 2.4656222608, 2.5800391630], [2.1352915140, 1.2328111304, 6.0071852168]])
plt.rcParams.update({"font.size": 8.5, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans"})


def spin(order, i, layer):
    sub = 1 if layer == 0 else -1
    return sub if order == "A" else sub * (-1) ** i


fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.75), gridspec_kw={"width_ratios": [0.8, 1, 1], "wspace": 0.06})

# (a) side view: two cells along x, two along c, projected along a2-ish (depth = y)
a = ax[0]
atoms = []
for i, k in itertools.product(range(-1, 3), range(0, 2)):
    t = i * A[0] + k * A[2]
    for l in (0, 1):
        atoms.append(("Mn", l, MN[l] + t))
        atoms.append(("Se", l, SEL[l] + t))
xlim = (-0.8, 2 * 4.27 + 1.2)
zlim = (-0.9, 2 * 6.854 + 0.4)
vis = [x for x in atoms if xlim[0] < x[2][0] < xlim[1] and zlim[0] < x[2][2] < zlim[1]]
for p, q in itertools.combinations(vis, 2):
    if {p[0], q[0]} == {"Mn", "Se"}:
        d = q[2] - p[2]
        if abs(d[1]) < 1e-3 or np.linalg.norm(d) < 2.72:
            if np.linalg.norm(d) < 2.72:
                a.plot([p[2][0], q[2][0]], [p[2][2], q[2][2]], color="#C9CCD2", lw=0.9, zorder=1)
for el, l, r in sorted(vis, key=lambda x: -x[2][1]):
    if el == "Se":
        a.scatter(r[0], r[2], s=26, color=SE, edgecolor="white", lw=0.4, zorder=2)
    else:
        up = l == 0
        a.scatter(r[0], r[2], s=62, marker="^" if up else "v", color=UP if up else DN, edgecolor="white", lw=0.4, zorder=3)
a.annotate("", xy=(xlim[1] - 0.4, 2.4 + 3.2), xytext=(xlim[1] - 0.4, 2.4), arrowprops=dict(arrowstyle="->", color=INK, lw=0.8))
a.text(xlim[1] - 0.15, 4.0, "$c$", color=INK, va="center", ha="left")
a.set_xlim(xlim[0], xlim[1] + 0.8)
a.set_ylim(*zlim)
a.set_aspect("equal")
a.axis("off")
a.set_title("(a)", loc="left", fontsize=8.5, color=INK)

# (b, c) view along c
for axx, order, title in ((ax[1], "A", "(b) A-type"), (ax[2], "S", r"(c) stripe AFMs$\alpha$")):
    for i, j in itertools.product(range(-4, 9), range(-1, 7)):
        t = i * A[0] + j * A[1]
        for l in (1, 0):
            r = MN[l] + t
            if not (-0.6 < r[0] < 13.5 and -0.3 < r[1] < 14.2):
                continue
            s = spin(order, i, l)
            axx.scatter(r[0], r[1], s=70 if l == 0 else 24, marker="^" if s > 0 else "v", color=UP if s > 0 else DN,
                        edgecolor="white", lw=0.4, zorder=3 if l == 0 else 2)
    n = 1 if order == "A" else 2
    o = MN[0] + 1 * A[0] + 1 * A[1]
    v = [o, o + n * A[0], o + n * A[0] + A[1], o + A[1], o]
    axx.plot([p[0] for p in v], [p[1] for p in v], color=INK, lw=0.8, ls=(0, (3, 2)), zorder=4)
    axx.set_aspect("equal")
    axx.axis("off")
    axx.set_title(title, loc="left", fontsize=8.5, color=INK)
fig.subplots_adjust(left=0.01, right=0.99, top=0.9, bottom=0.03)
fig.savefig(OUT + ".pdf")
fig.savefig(OUT + ".png", dpi=300)
print("wrote", os.path.abspath(OUT))
