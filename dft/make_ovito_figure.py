#!/usr/bin/env python3
"""Assemble the OVITO renders (figures/ovito/*.png, from render_ovito.sh) into one labelled figure.
Writes figures/fig_orders_ovito.{png,pdf}."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from PIL import Image, ImageChops  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = os.path.join(HERE, "..", "figures")
UP, DN, SE, INK = "#2E5FA3", "#C1442E", "#BDC2C9", "#1F2328"
PANELS = [("rs_AFMII", r"(a) $\alpha$-MnSe, AFM-II"), ("wz_Atype", "(b) w-MnSe, A-type"),
          ("wz_AFMsa", r"(c) w-MnSe, stripe AFMs$\alpha$")]
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans"})


def load(name):
    i = Image.open(os.path.join(FIG, "ovito", name + ".png")).convert("RGB")
    return np.asarray(i.crop(ImageChops.difference(i, Image.new("RGB", i.size, "white")).getbbox()))


imgs = [load(n) for n, _ in PANELS]
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.95), gridspec_kw={"width_ratios": [i.shape[1] / i.shape[0] for i in imgs], "wspace": 0.04})
for a, im, (_, t) in zip(ax, imgs, PANELS):
    a.imshow(im, interpolation="lanczos")
    a.set_axis_off()
    a.set_title(t, loc="left", fontsize=9, color=INK, pad=4)
h = [Line2D([], [], marker="o", ls="", color=UP, ms=7, label=r"Mn, spin $\uparrow$"),
     Line2D([], [], marker="o", ls="", color=DN, ms=7, label=r"Mn, spin $\downarrow$"),
     Line2D([], [], marker="o", ls="", color=SE, ms=5, label="Se")]
fig.legend(handles=h, loc="lower center", ncol=3, frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, -0.005))
fig.subplots_adjust(left=0.005, right=0.995, top=0.91, bottom=0.085)
for ext, dpi in (("png", 400), ("pdf", 400)):
    fig.savefig(os.path.join(FIG, "fig_orders_ovito." + ext), dpi=dpi, facecolor="white")
print("wrote", os.path.abspath(os.path.join(FIG, "fig_orders_ovito.png")))
