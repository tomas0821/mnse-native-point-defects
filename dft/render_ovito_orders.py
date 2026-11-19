#!/usr/bin/env python3
"""OVITO renderings of the magnetic orders: AFM-II rock-salt MnSe, and A-type and stripe wurtzite MnSe.

Geometries are the relaxed cells used in the calculations. Mn colour and arrow give the spin. The spins are
collinear and the calculations have no spin-orbit coupling, so the arrow direction relative to the lattice is
schematic: arrows are drawn along the page's vertical for readability, along c in wurtzite and [111] in rock salt.
Each panel carries axes: a1, a2, c for wurtzite; [111] and [-1-12] for rock salt.

usage: ~/.venvs/ovito/bin/python dft/render_ovito_orders.py [--preview] [scene ...]
scenes: wz_Atype wz_AFMsa wz_AFMsb rs_AFMII. Writes figures/ovito/<scene>.png
"""
import itertools
import os
import re
import sys

import numpy as np
from ovito.data import DataCollection, Particles, ParticleType, SimulationCell
from ovito.modifiers import CreateBondsModifier
from ovito.pipeline import Pipeline, StaticSource
from ovito.qt_compat import QtCore
from ovito.vis import BondsVis, CoordinateTripodOverlay, TachyonRenderer, VectorVis, Viewport

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures", "ovito")
UP = (0x2E / 255, 0x5F / 255, 0xA3 / 255)
DN = (0xC1 / 255, 0x44 / 255, 0x2E / 255)
SE = (0.74, 0.76, 0.79)
A = np.array([[4.270583028, 0, 0], [-2.135291514, 3.698433391, 0], [0, 0, 6.854295104]])
MN = np.array([[0.0, 2.4656222608, -0.0096777392], [2.1352915140, 1.2328111304, 3.4174697914]])
SEL = np.array([[0.0, 2.4656222608, 2.5800391630], [2.1352915140, 1.2328111304, 6.0071852168]])


def wz_spin(order, i, layer):
    sub = 1 if layer == 0 else -1
    if order == "Atype":
        return sub
    if order == "AFMsa":
        return sub * (-1) ** i
    return (-1) ** i          # AFMsb: same chain pattern in both layers


def wurtzite(order, rep=(4, 2, 1)):
    """Block of the lattice, with atoms on the far faces included so that the block looks closed."""
    pos, el, spin = [], [], []
    for i, j, k in itertools.product(range(-1, rep[0] + 1), range(-1, rep[1] + 1), range(-1, rep[2] + 1)):
        t = i * A[0] + j * A[1] + k * A[2]
        for l in (0, 1):
            for r, e in ((MN[l] + t, "Mn"), (SEL[l] + t, "Se")):
                f = np.linalg.solve(A.T, r)
                if np.all(f > -0.02) and f[0] < rep[0] + 0.02 and f[1] < rep[1] + 0.02 and f[2] < rep[2] + 0.02:
                    pos.append(r)
                    el.append(e)
                    spin.append(wz_spin(order, i, l) if e == "Mn" else 0)
    return np.array(pos), el, spin


def rocksalt():
    t = open(os.path.join(HERE, "prod", "d_rs_pristine_s120.in")).read()
    cell = np.array([[float(x) for x in l.split()] for l in
                     re.search(r"CELL_PARAMETERS angstrom\n((?:.*\n){3})", t).group(1).splitlines()])
    body = t[t.index("ATOMIC_POSITIONS"):t.index("K_POINTS")].splitlines()[1:]
    pos, el, spin = [], [], []
    for line in body:
        p = line.split()
        if len(p) < 4:
            continue
        r0 = np.array(list(map(float, p[1:4])))
        f0 = np.linalg.solve(cell.T, r0)
        for sh in itertools.product((0, 1), repeat=3):
            f = f0 + sh
            if np.all(f < 1.02) and np.all(f > -0.02):
                pos.append(cell.T @ f)
                el.append("Se" if p[0] == "Se" else "Mn")
                spin.append(0 if p[0] == "Se" else (1 if p[0] == "Mn1" else -1))
    return np.array(pos), el, spin, cell


def tripod(name):
    """Axes in neutral ink, so they cannot be read as spin colours. Wurtzite: the lattice vectors a1, a2, c.
    Rock salt: [111], the stacking direction of the ferromagnetic sheets (vertical here), and [-1-12] across the page."""
    ink = (0.16, 0.17, 0.19)
    t = CoordinateTripodOverlay(size=0.10, line_width=0.06, font_size=0.42, style=CoordinateTripodOverlay.Style.Flat,
                                alignment=QtCore.Qt.AlignmentFlag.AlignLeft | QtCore.Qt.AlignmentFlag.AlignBottom,
                                offset_x=0.0 if name.startswith("wz") else 0.05,
                                offset_y=0.16 if name.startswith("wz") else 0.05)
    if name.startswith("wz"):
        axes = [((1, 0, 0), "a<sub>1</sub>"), ((-0.5, 0.8660254, 0), "a<sub>2</sub>"), ((0, 0, 1), "c"), None]
    else:
        # Only the two directions that lie in the page: [111] up, [-1-12] across. The cubic axes all project obliquely
        # in this view and two of them overlap.
        ob = "<span style='text-decoration: overline'>1</span>"
        axes = [((1 / 3 ** 0.5,) * 3, "[111]"), ((-1 / 6 ** 0.5, -1 / 6 ** 0.5, 2 / 6 ** 0.5), f"[{ob}&#8202;{ob}&#8202;2]"), None, None]
    for i, ax in enumerate(axes, start=1):
        if ax is None:
            setattr(t, f"axis{i}_enabled", False)
            continue
        d, lab = ax
        setattr(t, f"axis{i}_enabled", True)
        setattr(t, f"axis{i}_dir", d)
        setattr(t, f"axis{i}_label", lab)
        setattr(t, f"axis{i}_color", ink)
    return t


def scene(name, preview):
    if name.startswith("wz"):
        order = name.split("_")[1]
        pos, el, spin = wurtzite(order)
        n = 1 if order == "Atype" else 2
        mcell = np.column_stack([n * A[0], A[1], A[2]])
        origin = MN[0] + A[0]
        axis = np.array([0, 0, 1.0])
        cutoff = 2.85
        d = A[1] / np.linalg.norm(A[1])         # look along the stripe chains, with a small tilt
        view = (d[0] - 0.28, d[1] - 0.10, -0.20)
        up = (0, 0, 1)
    else:
        pos, el, spin, cell = rocksalt()
        mcell = cell.T
        origin = np.zeros(3)
        axis = np.array([1.0, 1.0, 1.0]) / np.sqrt(3)
        cutoff = 3.0
        view = (0.80, -1.0, -0.12)
        up = (1, 1, 1)
    types = np.array([3 if e == "Se" else (1 if s > 0 else 2) for e, s in zip(el, spin)])
    data = DataCollection()
    c = SimulationCell(pbc=(False, False, False))
    c[:, :3] = mcell
    c[:, 3] = origin
    c.vis.line_width = 0.035
    c.vis.rendering_color = (0.12, 0.14, 0.16)
    c.vis.enabled = name.startswith("wz")
    data.objects.append(c)
    p = Particles()
    data.objects.append(p)
    p.create_property("Position", data=pos)
    pt = p.create_property("Particle Type", data=types)
    pt.types.append(ParticleType(id=1, name="MnUp", color=UP, radius=0.50))
    pt.types.append(ParticleType(id=2, name="MnDown", color=DN, radius=0.50))
    pt.types.append(ParticleType(id=3, name="Se", color=SE, radius=0.32))
    vec = np.array([s * axis for s in spin], dtype=float)
    sp = p.create_property("SpinVector", dtype=float, components=3, data=vec)
    sp.vis = VectorVis(alignment=VectorVis.Alignment.Center, scaling=2.9, width=0.12, flat_shading=False)
    p.create_property("Vector Color", data=np.array([np.array(UP) * 0.72 if s > 0 else np.array(DN) * 0.72 for s in spin]))
    pipe = Pipeline(source=StaticSource(data=data))
    bonds = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise, vis=BondsVis(width=0.11, color=(0.78, 0.79, 0.81),
                                coloring_mode=BondsVis.ColoringMode.Uniform))
    bonds.set_pairwise_cutoff("MnUp", "Se", cutoff)
    bonds.set_pairwise_cutoff("MnDown", "Se", cutoff)
    pipe.modifiers.append(bonds)
    pipe.add_to_scene()
    vp = Viewport(type=Viewport.Type.Ortho, camera_dir=view, camera_up=up)
    size = (900, 760) if preview else (2600, 2200)
    vp.zoom_all(size)
    vp.fov *= 1.06                       # a little room in the corner for the axes
    vp.overlays.append(tripod(name))
    r = TachyonRenderer(ambient_occlusion=True, ambient_occlusion_brightness=0.75, direct_light_intensity=0.85,
                        shadows=True, antialiasing_samples=4 if preview else 24)
    f = os.path.join(OUT, name + ".png")
    vp.render_image(size=size, filename=f, background=(1, 1, 1), renderer=r)
    pipe.remove_from_scene()
    print("wrote", os.path.abspath(f), "atoms", len(pos))


def main():
    a = sys.argv[1:]
    preview = "--preview" in a
    names = [x for x in a if not x.startswith("--")] or ["wz_Atype", "wz_AFMsa", "wz_AFMsb", "rs_AFMII"]
    os.makedirs(OUT, exist_ok=True)
    for n in names:
        scene(n, preview)


if __name__ == "__main__":
    main()
