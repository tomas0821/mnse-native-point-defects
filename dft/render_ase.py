#!/usr/bin/env python3
"""
ASE-rendered structure figures, spin-resolved.

Two backends:
  --mpl   built-in matplotlib renderer. No extra dependencies. Proper depth
          sorting and sphere overlap, which the hand-rolled scatter plots in
          make_structure_figures.py do not give.
  --pov   writes POV-Ray scene files for ray-traced output (shadows, specular
          highlights, true publication quality). Needs the `povray` binary.

Colours are the CVD-validated palette used everywhere else in this project
(normal-vision dE_OKLab >= 15, worst-case dichromat >= 8). Spin is encoded by
colour AND by the Mn1/Mn2 species split, so it survives greyscale.

    python3 dft/render_ase.py --mpl  figures/structures/*.cif
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.visualize.plot import plot_atoms

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from qe2vesta import parse_qe          # single source of truth for geometry+moments


def load(path):
    """Read a QE output directly. CIF CANNOT carry magnetic moments -- going via
    CIF silently drops the spins and every Mn renders the same colour."""
    if path.endswith(".out"):
        cell, labels, pos, mag = parse_qe(path)
        syms = ["Mn" if l.startswith("Mn") else "Se" for l in labels]
        at = Atoms(symbols=syms, positions=pos, cell=cell, pbc=True)
        at.set_initial_magnetic_moments(mag)
        return at
    at = read(path)
    if abs(at.get_initial_magnetic_moments()).max() == 0:
        print(f"  !! {os.path.basename(path)}: no magnetic moments in this file "
              f"(CIF cannot store them) -- pass the QE .out instead")
    return at

UP, DN, SE = "#2E5FA3", "#C1442E", "#A8AEB8"
INK, MUTED = "#1F2328", "#6B7280"


def _rgb(h):
    """POV-Ray needs rgb floats; hex strings are passed through verbatim by ASE
    and produce `color #2E5FA3`, which POV-Ray cannot parse."""
    h = h.lstrip("#")
    return tuple(int(h[i:i + 2], 16) / 255 for i in (0, 2, 4))


def colour_and_radii(at):
    """Mn split by the sign of its magnetic moment; Se recessive and smaller."""
    mag = at.get_initial_magnetic_moments()
    cols, rad = [], []
    for s, m in zip(at.get_chemical_symbols(), mag):
        if s == "Mn":
            cols.append(_rgb(UP) if m >= 0 else _rgb(DN))
            rad.append(0.62)
        else:
            cols.append(_rgb(SE))
            rad.append(0.42)
    return cols, np.array(rad)


def render_mpl(paths, out, rotation="-72x,18y,0z"):
    n = len(paths)
    fig, axes = plt.subplots(1, n, figsize=(5.2 * n, 5.4), dpi=300)
    if n == 1:
        axes = [axes]
    for ax, p in zip(axes, paths):
        at = load(p)
        cols, rad = colour_and_radii(at)
        plot_atoms(at, ax, radii=rad, colors=cols, rotation=rotation,
                   show_unit_cell=2)
        tag = os.path.basename(p).replace(".cif","").replace(".out","")
        mag = at.get_initial_magnetic_moments()
        up = int((mag > 0.2).sum()); dn = int((mag < -0.2).sum())
        ax.set_title(f"{tag}\n{len(at)} atoms · Mn↑{up} Mn↓{dn}",
                     fontsize=10, color=INK, weight="bold")
        ax.set_axis_off()
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], marker="o", ls="", ms=10, mfc=UP, mec="white", label="Mn  spin ↑"),
                        Line2D([], [], marker="o", ls="", ms=10, mfc=DN, mec="white", label="Mn  spin ↓"),
                        Line2D([], [], marker="o", ls="", ms=8, mfc=SE, mec="white", label="Se")],
               loc="lower center", ncol=3, frameon=False, fontsize=10)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(out, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  wrote {out} ({os.path.getsize(out)//1024} kB)")


def render_pov(paths, outdir):
    os.makedirs(outdir, exist_ok=True)
    for p in paths:
        at = load(p)
        cols, rad = colour_and_radii(at)
        tag = os.path.basename(p).replace(".cif","").replace(".out","")
        dst = os.path.join(outdir, tag + ".pov")
        write(dst, at, format="pov", radii=rad, colors=cols,
              rotation="-72x,18y,0z",
              povray_settings=dict(canvas_width=1800, transparent=False,
                                   camera_type="orthographic",
                                   background="White", display=False))
        print(f"  wrote {dst}  ->  render with:  povray +I{tag}.ini")


if __name__ == "__main__":
    args = sys.argv[1:]
    mode = "--mpl"
    if args and args[0].startswith("--"):
        mode, args = args[0], args[1:]
    if not args:
        sys.exit("usage: render_ase.py [--mpl|--pov] <structure.cif> ...")
    if mode == "--mpl":
        render_mpl(args, "figures/fig5_ase_structures.png")
    else:
        render_pov(args, "figures/pov")
