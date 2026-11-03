#!/usr/bin/env python3
"""Re-derive the Stage 4 local-structure / local-moment numbers from the per-site .dat files
(INIT/FINAL Cartesian positions in Angstrom, MAGN sphere-integrated site moments), locating the
defect by comparison with the pristine cell rather than by re-running the builders.
usage: local_analysis.py   (run inside dft/analysis)"""
import re
import numpy as np


def load(fn):
    sec, d = None, {"CELL": [], "INIT": [], "FINAL": [], "MAGN": []}
    for line in open(fn):
        if line.startswith("#"):
            sec = line[1:].split()[0]
            continue
        if sec in ("CELL",):
            d["CELL"].append([float(x) for x in line.split()])
        elif sec in ("INIT", "FINAL"):
            p = line.split()
            d[sec].append((re.sub(r"\d", "", p[0]), np.array([float(x) for x in p[1:4]])))
        elif sec == "MAGN":
            m = re.search(r"magn=\s*(-?\d+\.\d+)", line)
            if m:
                d["MAGN"].append(float(m.group(1)))
    d["CELL"] = np.array(d["CELL"])
    return d


def qe_positions(fn):
    txt = open(fn).read()
    blk = txt[txt.index("ATOMIC_POSITIONS"):].splitlines()[1:]
    out = []
    for line in blk:
        p = line.split()
        if len(p) < 4 or not re.match(r"[A-Z]", p[0]):
            break
        out.append((re.sub(r"\d", "", p[0]), np.array([float(x) for x in p[1:4]])))
    return out


def mic(v, cell):
    f = np.linalg.solve(cell.T, v)
    f -= np.round(f)
    return cell.T @ f


def analyse(pol, defect):
    P = load(f"d_{pol}_prist.dat")
    if not P["INIT"]:                     # pristine .dat holds moments only -> take positions
        P["INIT"] = qe_positions(f"../prod/d_{pol}_pristine_s120.in")   # from the 120 Ry input
    D = load(f"d_{pol}_{defect}.dat")
    cell = D["CELL"]
    # match every pristine site to the nearest defect-cell INIT atom
    dpos = [p for _, p in D["INIT"]]
    used, site = set(), None
    for i, (sp, p) in enumerate(P["INIT"]):
        dist = [np.linalg.norm(mic(q - p, cell)) for q in dpos]
        j = int(np.argmin(dist))
        if dist[j] > 0.3:                      # pristine site with no atom -> vacancy
            site = (i, sp, p, None)
        elif D["INIT"][j][0] != sp:            # species changed -> antisite
            site = (i, sp, p, j)
        used.add(j)
    assert site is not None, "defect site not found"
    i0, sp0, p0, j0 = site
    pm = np.array(P["MAGN"])
    ref = np.median(np.abs(pm[[k for k, (s, _) in enumerate(P["INIT"]) if s == "Mn"]]))
    rows = []
    for j, ((sp, pi), (_, pf), m) in enumerate(zip(D["INIT"], D["FINAL"], D["MAGN"])):
        r = mic(pi - p0, cell)
        R = np.linalg.norm(r)
        disp = mic(pf - pi, cell)
        rows.append((j, sp, R, float(disp @ r / R) if R > 0.1 else 0.0, np.linalg.norm(disp), m))
    shell = sorted(x[2] for x in rows if x[2] > 0.1)[0]
    first = [x for x in rows if x[2] > 0.1 and abs(x[2] - shell) < 0.15]
    mn = [x for x in rows if x[1] == "Mn" and x[0] != j0]
    dm = [(abs(abs(x[5]) - ref), x[2]) for x in mn]
    big = max(dm)
    out = dict(pol=pol, defect=defect, shell_n=len(first), shell_sp=first[0][1], shell_R=shell,
               radial=np.mean([x[3] for x in first]), maxdisp=max(x[4] for x in rows),
               ref=ref, max_dm=big[0], max_dm_R=big[1])
    if j0 is not None and D["INIT"][j0][0] == "Mn":          # Mn antisite: its own moment
        out["antisite_m"] = D["MAGN"][j0]
        out["shell_moments"] = [round(x[5], 2) for x in first]
    return out


if __name__ == "__main__":
    for pol in ("rs", "wz"):
        for d in ("V_Mn", "V_Se", "Se_Mn", "Mn_Se"):
            o = analyse(pol, d)
            line = (f"{pol} {d:6s} 1st shell {o['shell_n']} {o['shell_sp']:2s} @ {o['shell_R']:.3f} A"
                    f" | radial {o['radial']:+.3f} A | max|disp| {o['maxdisp']:.3f} A"
                    f" | bulk |m| {o['ref']:.3f} | max|dm| {o['max_dm']:.3f} at R={o['max_dm_R']:.2f} A")
            if "antisite_m" in o:
                line += f" | antisite Mn m = {o['antisite_m']:+.3f}, 1st-shell m = {o['shell_moments']}"
            print(line)
