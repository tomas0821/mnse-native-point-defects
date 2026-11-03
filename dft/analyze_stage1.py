#!/usr/bin/env python3
"""
Turn the Stage 1 outputs into decisions.

Run from the project root once the batch has finished and the .out files are back:

    python3 dft/analyze_stage1.py dft/stage1

DECISION RULES ENCODED HERE (see 00_Master_Notebook.md "How the values get decided"):

  Cutoff   : smallest ecutwfc whose energy/atom is within TOL_E of the 90 Ry
             reference AND whose pressure is within TOL_P. Pressure matters
             because Batch C is vc-relax -- an under-converged cutoff produces
             Pulay stress and a wrong lattice constant.
  k-mesh   : smallest mesh within TOL_E/atom of the densest mesh.
  Orderings: report dE per formula unit between magnetic orderings. This is the
             quantity the paper actually uses, and it converges far faster than
             the absolute energy -- so it is also reported at every cutoff as a
             cross-check that the *difference* is converged, not just E.

Thresholds are deliberately tighter than the physics needs, because the errors
propagate into formation energies (target ~10 meV) and into dE between orderings
(target ~1 meV/f.u.).
"""
import os
import re
import sys
import glob

RY = 13.605693122994          # eV
TOL_E = 1.0                   # meV/atom
TOL_P = 0.5                   # kbar


def parse(path):
    txt = open(path, errors="ignore").read()
    if "JOB DONE" not in txt:
        return None
    d = {}
    e = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", txt, re.M)
    if not e:
        return None
    d["E_Ry"] = float(e[-1])
    n = re.search(r"number of atoms/cell\s+=\s+(\d+)", txt)
    d["nat"] = int(n.group(1)) if n else None
    p = re.findall(r"P=\s*(-?\d+\.\d+)", txt)
    d["P_kbar"] = float(p[-1]) if p else None
    f = re.findall(r"Total force =\s+(\d+\.\d+)", txt)
    d["F"] = float(f[-1]) if f else None
    am = re.findall(r"absolute magnetization\s+=\s+(-?\d+\.\d+)", txt)
    d["absmag"] = float(am[-1]) if am else None
    tm = re.findall(r"total magnetization\s+=\s+(-?\d+\.\d+)", txt)
    d["totmag"] = float(tm[-1]) if tm else None
    # final cell for vc-relax
    cp = re.findall(r"CELL_PARAMETERS \(angstrom\)\n(.*?)\n\n", txt, re.S)
    if cp:
        d["cell"] = [[float(x) for x in row.split()]
                     for row in cp[-1].strip().splitlines()]
    d["nsteps"] = len(e)
    return d


def emev_per_atom(d):
    return d["E_Ry"] * RY * 1000.0 / d["nat"]


def table(rows, headers):
    w = [max(len(str(h)), max((len(str(r[i])) for r in rows), default=0))
         for i, h in enumerate(headers)]
    print("  " + "  ".join(str(h).ljust(w[i]) for i, h in enumerate(headers)))
    print("  " + "  ".join("-" * w[i] for i in range(len(headers))))
    for r in rows:
        print("  " + "  ".join(str(c).ljust(w[i]) for i, c in enumerate(r)))


def scan(d_, tag_re, label, xkey):
    got = []
    for p in sorted(glob.glob(os.path.join(d_, "*.out"))):
        m = re.match(tag_re, os.path.basename(p))
        if not m:
            continue
        r = parse(p)
        if r:
            got.append((int(m.group(1)), r))
    if not got:
        print(f"\n### {label}: no completed runs yet\n")
        return None
    got.sort()
    ref = got[-1][1]
    rows = []
    pick = None
    for x, r in got:
        de = abs(emev_per_atom(r) - emev_per_atom(ref))
        dp = abs((r["P_kbar"] or 0) - (ref["P_kbar"] or 0))
        ok = de < TOL_E and dp < TOL_P
        if ok and pick is None:
            pick = x
        rows.append([x, f"{r['E_Ry']:.8f}", f"{de:8.3f}",
                     f"{r['P_kbar']:.2f}" if r["P_kbar"] is not None else "-",
                     f"{dp:6.2f}", f"{r['absmag']:.3f}" if r["absmag"] else "-",
                     "OK" if ok else ""])
    print(f"\n### {label}   (reference = densest/highest: {got[-1][0]})")
    table(rows, [xkey, "E (Ry)", "dE meV/at", "P kbar", "dP", "|m|", ""])

    # ---- PLATEAU TEST -------------------------------------------------------
    # Measuring every point against the densest one is NOT sufficient: the
    # reference is trivially within 0 of itself, so it always looks "OK" and the
    # scan reports convergence even when it never plateaued. Require the LAST TWO
    # points to agree before trusting any pick.
    if len(got) < 2:
        print("  --> only one point; cannot assess convergence.")
        return None
    tail_de = abs(emev_per_atom(got[-1][1]) - emev_per_atom(got[-2][1]))
    tail_dp = abs((got[-1][1]["P_kbar"] or 0) - (got[-2][1]["P_kbar"] or 0))
    print(f"  plateau test: last two points ({got[-2][0]} -> {got[-1][0]}) differ by "
          f"{tail_de:.3f} meV/atom, {tail_dp:.2f} kbar")
    if tail_de >= TOL_E or tail_dp >= TOL_P:
        print(f"  --> !! NOT CONVERGED. The scan has not plateaued -- the largest value "
              f"tested ({got[-1][0]}) is still moving.")
        print(f"     Any 'OK' above is measured against a moving reference and is "
              f"MEANINGLESS. EXTEND THE SCAN beyond {got[-1][0]}.")
        return None
    if pick is not None:
        print(f"  --> CONVERGED CHOICE: {xkey} = {pick} "
              f"(first within {TOL_E} meV/atom and {TOL_P} kbar of a plateaued reference)")
    return pick


def orderings(d_):
    print("\n### Magnetic ordering energetics")
    for sysname, fu_atoms in (("rs", 2), ("wz", 2)):
        got = {}
        for p in glob.glob(os.path.join(d_, f"ordB_{sysname}_*.out")):
            tag = re.search(rf"ordB_{sysname}_(\w+)\.out", os.path.basename(p)).group(1)
            r = parse(p)
            if r:
                got[tag] = r
        if not got:
            continue
        nfu = list(got.values())[0]["nat"] / fu_atoms
        base = min(got, key=lambda k: got[k]["E_Ry"])
        rows = []
        for tag, r in sorted(got.items(), key=lambda kv: kv[1]["E_Ry"]):
            dE = (r["E_Ry"] - got[base]["E_Ry"]) * RY * 1000.0 / nfu
            rows.append([tag, f"{r['E_Ry']:.8f}", f"{dE:9.2f}",
                         f"{r['totmag']:.3f}" if r["totmag"] is not None else "-",
                         f"{r['absmag']:.3f}" if r["absmag"] else "-",
                         "<-- ground state" if tag == base else ""])
        print(f"\n  {sysname.upper()}   ({int(nfu)} formula units/cell)")
        table(rows, ["ordering", "E (Ry)", "dE meV/f.u.", "M_tot", "|m|", ""])
        if sysname == "rs":
            if base != "AFMII":
                print("  !! WARNING: ground state is not AFM-II. Expected type-II for "
                      "alpha-MnSe (fcc type-II antiferromagnet, like MnO/NiO).")
                print("     Check U, the starting magnetisations, and that the SCF "
                      "did not collapse to a different state than it started in.")
            else:
                print("  OK: AFM-II is the ground state, as expected.")


def relax(d_):
    print("\n### vc-relax results")
    for p in sorted(glob.glob(os.path.join(d_, "relC_*.out"))):
        r = parse(p)
        name = os.path.basename(p).replace(".out", "")
        if not r:
            print(f"  {name}: not finished")
            continue
        print(f"\n  {name}  ({r['nsteps']} SCF cycles)")
        if "cell" in r:
            import math
            c = r["cell"]
            L = [math.dist([0, 0, 0], v) for v in c]
            ang = []
            for i, j in ((1, 2), (0, 2), (0, 1)):
                cs = sum(c[i][k] * c[j][k] for k in range(3)) / (L[i] * L[j])
                ang.append(math.degrees(math.acos(max(-1, min(1, cs)))))
            print(f"    |a1|,|a2|,|a3| = {L[0]:.4f} {L[1]:.4f} {L[2]:.4f} A")
            print(f"    angles         = {ang[0]:.3f} {ang[1]:.3f} {ang[2]:.3f} deg")
            if name.startswith("relC_rs"):
                dev = max(abs(a - 90.0) for a in ang)
                n = int(round(L[0] / 5.4616))
                print(f"    implied cubic a = {L[0]/max(n,1):.4f} A "
                      f"(experiment 5.4616 A, Pradhan 2025)")
                if dev > 0.05:
                    print(f"    ** RHOMBOHEDRAL DISTORTION: angles deviate by {dev:.3f} deg "
                          f"from 90. This is the magnetostructural coupling -- it means the "
                          f"defect supercells must NOT assume an undistorted cubic cell.")
                else:
                    print("    no significant distortion (angles within 0.05 deg of 90)")
        print(f"    |m| = {r['absmag']}  M_tot = {r['totmag']}")


if __name__ == "__main__":
    d = sys.argv[1] if len(sys.argv) > 1 else "dft/stage1"
    done = len(glob.glob(os.path.join(d, "*.done")))
    outs = len(glob.glob(os.path.join(d, "*.out")))
    print(f"Stage 1 analysis: {done} .done, {outs} .out in {d}")
    print(f"Thresholds: {TOL_E} meV/atom, {TOL_P} kbar")
    scan(d, r"convA_rs_ecut(\d+)\.out", "alpha-MnSe cutoff scan", "ecutwfc")
    scan(d, r"convA_wz_ecut(\d+)\.out", "w-MnSe cutoff scan", "ecutwfc")
    scan(d, r"convA_rs_k(\d+)\.out", "alpha-MnSe k-mesh scan", "k")
    scan(d, r"convA_wz_k(\d+)\.out", "w-MnSe k-mesh scan", "k")
    orderings(d)
    relax(d)
    print("\nNOTE: U is NOT decided by this script -- it is not a convergence "
          "parameter. See the notebook section 'How the values get decided'.")
