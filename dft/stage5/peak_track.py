#!/usr/bin/env python3
"""Independent cross-check of unfold_analyse.py: track the DOMINANT unfolded peak near the
pristine band instead of stacking cumulative weight.

For each spin and K: broaden the unfolded spectrum A(E) = sum_n w_n g(E - E_n; sigma), and within
+-dE of the pristine primitive band energy (shifted by the Se-4s alignment) take the weight-averaged
energy of the states around the maximum of A(E). This needs no global band counting, so it is
insensitive to the unequal spectral-weight loss between spins that a net moment causes.

usage: peak_track.py WEIGHTS.npz BANDS.out [--band 27] [--dE 0.20] [--sigma 0.03] [--dump OUT.npz]
Pristine reference: w_pristine.npz + wzsc_bands.out (states with w > 0.5, sorted, band index).
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from splitting import parse  # noqa: E402

REF_SE4S = -6.0521


def load(wf, bf):
    w = np.load(wf)["w"]
    b = parse(bf)
    return w, np.array([b["up"][1], b["dn"][1]])


def se4s(E, w):
    m = (E > -9.0) & (E < -3.0)
    return float((E[m] * w[m]).sum() / w[m].sum())


def peak(E, w, e0, dE, sig, hmin=0.35):
    """Nearest local maximum of the broadened spectrum to e0, of height at least hmin (in units
    where one fully coherent band has height 1). Nearest, not tallest: a taller neighbouring band
    inside the window must not capture the track."""
    grid = np.arange(e0 - dE, e0 + dE + 1e-9, 0.002)
    m = np.abs(E - e0) < dE + 4 * sig
    if not m.any():
        return np.nan, 0.0
    A = (w[m][:, None] * np.exp(-0.5 * ((grid[None, :] - E[m][:, None]) / sig) ** 2)).sum(0)
    i = np.where((A[1:-1] >= A[:-2]) & (A[1:-1] > A[2:]) & (A[1:-1] >= hmin))[0] + 1
    if len(i) == 0:
        return np.nan, 0.0
    ep = grid[i[np.argmin(np.abs(grid[i] - e0))]]
    n = np.abs(E - ep) < 2.5 * sig
    return (float((E[n] * w[n]).sum() / w[n].sum()), float(w[n].sum())) if w[n].sum() > 0 else (ep, 0.0)


def compute(wf, bf, jb=27, dE=0.20, sig=0.03):
    """Track primitive band jb along the path; returns S (meV, 41 k), odd, even (21 t), wgt, shift."""
    here = os.path.dirname(os.path.abspath(__file__))
    wp, Ep = load(os.path.join(here, "w_pristine.npz"), os.path.join(here, "wzsc_bands.out"))
    w, E = load(wf, bf)
    shift = se4s(E, w) - REF_SE4S
    nk = E.shape[1]
    ref = np.zeros((2, nk))
    for s in range(2):
        for k in range(nk):
            sel = np.where((wp[s, k] > 0.5) & (Ep[s, k] < 6.82) & (Ep[s, k] > -2.5))[0]
            ref[s, k] = Ep[s, k][sel[jb - 21]] if len(sel) >= jb - 20 else np.nan
    eff = np.zeros((2, nk))
    wgt = np.zeros((2, nk))
    for s in range(2):
        for k in range(nk):
            eff[s, k], wgt[s, k] = peak(E[s, k], w[s, k], ref[s, k] + shift, dE, sig)
    S = (eff[0] - eff[1]) * 1e3
    t = np.arange(0, 21)
    return dict(S=S, odd=(S[20 - t] - S[20 + t]) / 2, even=(S[20 - t] + S[20 + t]) / 2, wgt=wgt,
                eff=eff, ref=ref, shift=shift)


def main():
    a = sys.argv[1:]
    opt = lambda k, d: type(d)(a[a.index(k) + 1]) if k in a else d
    jb, dE, sig = opt("--band", 27), opt("--dE", 0.20), opt("--sigma", 0.03)
    r = compute(a[0], a[1], jb, dE, sig)
    S, odd, even, wgt, shift = r["S"], r["odd"], r["even"], r["wgt"], r["shift"]
    print(f"alignment shift {shift*1e3:+.1f} meV; band {jb}; window +-{dE} eV, sigma {sig}")
    print(f"  peak weight along path: min {np.nanmin(wgt):.2f}, median {np.nanmedian(wgt):.2f}")
    print(f"  max |odd| = {np.nanmax(np.abs(odd)):.1f} meV at t={int(np.nanargmax(np.abs(odd)))} | "
          f"max |even| = {np.nanmax(np.abs(even)):.1f} meV | even at t=13: {even[13]:+.1f} | "
          f"s(L1,G,L2) = {S[0]:+.1f}, {S[20]:+.1f}, {S[40]:+.1f}")
    print(f"  at the peak t=13: weights up {wgt[0,7]:.2f}/{wgt[0,33]:.2f}, dn {wgt[1,7]:.2f}/{wgt[1,33]:.2f}; "
          f"tracked {np.isfinite(S).sum()}/{len(S)} k-points")
    if "--dump" in a:
        np.savez(a[a.index("--dump") + 1], **r)


if __name__ == "__main__":
    main()
