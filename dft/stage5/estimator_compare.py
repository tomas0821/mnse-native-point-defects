#!/usr/bin/env python3
"""Compare ways of reading band 27 off the unfolded spectrum, because the tracked-peak curves spike where a band
fragments. Estimators: tracked peak (peak_track.py); centroid in a hard window; centroid in a Gaussian window centred
on the pristine energy; the same with the centre iterated to self-consistency.
Validation on the pristine cell, and a bias test: move the reference by +d (up) and -d (down) and measure how much of
that artificial 2d leaks into the result. A centroid with a narrow window is pulled towards its reference.
usage: estimator_compare.py"""
import sys  # noqa: F401
import numpy as np
import peak_track as PT

REL = np.load("s_relaxed_u.npz")["S"]
OK = np.abs(np.arange(41) - 20) >= 4
T = np.arange(21)
odd = lambda S: (S[20 - T] - S[20 + T]) / 2


def soft(w, E, ref, shift, tau, it=0):
    c = np.zeros((2, 41))
    for s in range(2):
        for k in range(41):
            e0 = ref[s, k] + shift
            for _ in range(it + 1):
                g = w[s, k] * np.exp(-0.5 * ((E[s, k] - e0) / tau) ** 2)
                e0 = (E[s, k] * g).sum() / g.sum()
            c[s, k] = e0
    return (c[0] - c[1]) * 1e3


def hard(w, E, ref, shift, dE):
    c = np.zeros((2, 41))
    for s in range(2):
        for k in range(41):
            m = np.abs(E[s, k] - (ref[s, k] + shift)) < dE
            c[s, k] = (E[s, k][m] * w[s, k][m]).sum() / w[s, k][m].sum() if w[s, k][m].sum() > 0 else np.nan
    return (c[0] - c[1]) * 1e3


def stats(S):
    d = S[OK] - REL[OK]
    return f"peak {abs(odd(S)[13]):4.0f}  rms {np.sqrt(np.nanmean(d ** 2)):3.0f}  max {np.nanmax(np.abs(d)):3.0f}"


def tex_rows():
    """Rows for Table S14: peak, rms and largest deviation (meV) per defect and estimator."""
    out = []
    TEX = {"V_Mn": r"$V_{\mathrm{Mn}}$", "V_Se": r"$V_{\mathrm{Se}}$", "Se_Mn": r"$\mathrm{Se}_{\mathrm{Mn}}$",
           "Mn_Se": r"$\mathrm{Mn}_{\mathrm{Se}}$"}
    for d in TEX:
        wf, bf = f"w_{d}.npz", f"wzsc_{d}_bands.out"
        w, E = PT.load(wf, bf)
        r = PT.compute(wf, bf, 27, 0.20, 0.03)
        base = soft(w, E, r["ref"], r["shift"], 0.06)
        ref2 = r["ref"].copy()
        ref2[0] += 0.05
        ref2[1] -= 0.05
        leak = 100 * np.nanmean((soft(w, E, ref2, r["shift"], 0.06)[OK] - base[OK]) / 100.0)
        for name, S in (("Tracked peak, $\\sigma = 0.03$~eV", r["S"]),
                        ("Hard window, $\\pm0.10$~eV", hard(w, E, r["ref"], r["shift"], 0.10)),
                        ("Hard window, $\\pm0.15$~eV", hard(w, E, r["ref"], r["shift"], 0.15)),
                        ("Hard window, $\\pm0.17$~eV", hard(w, E, r["ref"], r["shift"], 0.17)),
                        ("Gaussian window, $\\tau = 0.06$~eV", base),
                        ("Gaussian window, self-consistent center", soft(w, E, r["ref"], r["shift"], 0.06, 5))):
            q = S[OK] - REL[OK]
            lk = f"{leak:.0f}\\%" if "tau" in name else "--"
            out.append(f"{TEX[d]} | {name} | {abs(odd(S)[13]):.0f} | {np.sqrt(np.nanmean(q ** 2)):.0f} | {np.nanmax(np.abs(q)):.0f} | {lk}")
    return out


def main():
    wp, Ep = PT.load("w_pristine.npz", "wzsc_bands.out")
    rp = PT.compute("w_pristine.npz", "wzsc_bands.out", 27, 0.20, 0.03)
    print("pristine, rms from the exact curve (meV):")
    for dE in (0.15, 0.18, 0.19, 0.20):
        print(f"  hard +-{dE:.2f} eV: {np.sqrt(np.nanmean((hard(wp, Ep, rp['ref'], rp['shift'], dE)[OK] - rp['S'][OK]) ** 2)):.1f}")
    for tau in (0.04, 0.06, 0.08):
        print(f"  soft tau {tau:.2f} eV: {np.sqrt(np.nanmean((soft(wp, Ep, rp['ref'], rp['shift'], tau)[OK] - rp['S'][OK]) ** 2)):.2f}")
    for row in tex_rows():
        print(row)


if __name__ == "__main__":
    main()
