#!/usr/bin/env python3
"""Spin splitting from unfolded supercell spectra (Stage 5).

NOT used for the defect cells: stacking cumulative weight misassigns the two spins whenever a
defect removes spectral weight unequally between them (V_Mn: 10.72 vs 10.99), giving a spurious
uniform splitting. Use peak_track.py. Kept because it validates the pristine weights.

usage: unfold_analyse.py WEIGHTS.npz BANDS.out [--lo EV] [--hi EV] [--nb 11] [--band 7]
                         [--pristine] [--dump FILE]

Effective bands: within the valence manifold window [lo, hi] (default: the 11 primitive
bands 21-31 that sit above the Se-4s gap), each spin's unfolded spectral weight at each K is
sorted by energy, its cumulative sum is rescaled to nb, and effective band j is the
weight-averaged energy of the j-th unit of weight. For a pristine supercell (w = 0 or 1) this
reproduces the primitive bands exactly. Splitting s_j(K) = E_up_j - E_dn_j.

Path is L1 -(20)- Gamma -(20)- L2, 41 points, Gamma at index 20. Arms are compared at equal
distance t from Gamma: odd(t) = [s(20-t) - s(20+t)]/2 is the altermagnetic (sign-reversing)
part, even(t) = [s(20-t) + s(20+t)]/2 the uniform (ferrimagnetic-like) part. Pristine: even = 0.
"""
import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from splitting import parse  # noqa: E402


def effective_bands(E, w, lo, hi, nb):
    sel = (E >= lo) & (E <= hi)
    e, ww = E[sel], w[sel]
    o = np.argsort(e)
    e, ww = e[o], ww[o]
    tot = ww.sum()
    if tot <= 0:
        return np.full(nb, np.nan), 0.0
    c = np.concatenate([[0.0], np.cumsum(ww)]) * nb / tot
    out = np.empty(nb)
    for j in range(nb):
        a, b = j, j + 1
        ov = np.clip(np.minimum(c[1:], b) - np.maximum(c[:-1], a), 0, None)
        out[j] = (ov * e).sum() / ov.sum()
    return out, tot


def main():
    a = sys.argv[1:]
    d = np.load(a[0])
    w = d["w"]
    bands = parse(a[1])
    (ku, eu), (kd, ed) = bands["up"], bands["dn"]
    assert ku == kd
    E = np.array([eu, ed])                       # [spin, k, band] eV
    nk = E.shape[1]
    assert w.shape == E.shape, f"weights {w.shape} vs eigenvalues {E.shape}"
    xk_out = np.array(ku)
    # wfc files store xk in bohr^-1; the text output in 2pi/alat. Same order => constant ratio.
    nz = np.abs(xk_out) > 1e-3
    r = d["xk"][nz] / xk_out[nz]
    assert np.allclose(r, r.mean(), rtol=2e-3), "k-point order mismatch"
    assert np.all(np.abs(d["xk"][~nz]) < 1e-3 * r.mean() * 10), "k-point order mismatch (zeros)"
    opt = lambda k, dflt: type(dflt)(a[a.index(k) + 1]) if k in a else dflt
    nb, jb = opt("--nb", 11), opt("--band", 7) - 1
    # Window edges are fixed on the PRISTINE 72-atom spectrum (Se-4s gap -5.82..0.85 eV,
    # main gap 5.96..7.67 eV) and shifted by the potential alignment of this cell: the
    # weight-averaged energy of the Se-4s band (unfolded weights suppress localised defect
    # states). PRISTINE reference value REF_SE4S is measured by --pristine runs.
    REF_SE4S = opt("--ref", -6.0521)   # measured on pristine 72-atom cell, 2026-11-03
    se = (E > -9.0) & (E < -3.0)
    se4s = float((E[se] * w[se]).sum() / w[se].sum())
    shift = se4s - REF_SE4S if "--pristine" not in a else 0.0
    lo = opt("--lo", -2.5) + shift
    hi = opt("--hi", 6.82) + shift
    print(f"Se-4s weighted centroid = {se4s:.4f} eV; alignment shift = {shift*1e3:+.1f} meV")

    if "--pristine" in a:
        cnt = ((w[:, :, 414:558] > 0.5).sum(axis=2))
        frac01 = np.mean((w < 0.02) | (w > 0.98))
        print(f"[pristine] states with w>0.5 in window 415-558 per (spin,K): "
              f"min {cnt.min()} max {cnt.max()} (expect 8)")
        print(f"[pristine] fraction of all weights within 0.02 of 0 or 1: {frac01:.4f}")
        print(f"[pristine] use --ref {se4s:.4f} for defect cells")

    eff = np.zeros((2, nk, nb))
    tot = np.zeros((2, nk))
    for s in range(2):
        for k in range(nk):
            eff[s, k], tot[s, k] = effective_bands(E[s, k], w[s, k], lo, hi, nb)
    S = (eff[0] - eff[1]) * 1e3                  # meV, [k, band]
    sj = S[:, jb]
    t = np.arange(0, 21)
    odd = (sj[20 - t] - sj[20 + t]) / 2
    even = (sj[20 - t] + sj[20 + t]) / 2
    kp = int(np.argmax(np.abs(sj)))
    print(f"window [{lo:.2f}, {hi:.2f}] eV, nb={nb}; total weight per (spin,K): "
          f"{tot.min():.3f}-{tot.max():.3f}")
    print(f"effective band {jb+1} (primitive {jb+21}):")
    print(f"  peak |s| = {abs(sj[kp]):.1f} meV at K index {kp}")
    print(f"  arm1 / arm2 signed peak = {sj[:21][np.argmax(abs(sj[:21]))]:+.1f} / "
          f"{sj[20:][np.argmax(abs(sj[20:]))]:+.1f} meV")
    print(f"  s at L1 / Gamma / L2 = {sj[0]:+.1f} / {sj[20]:+.1f} / {sj[40]:+.1f} meV")
    print(f"  max |odd| (altermagnetic) = {np.max(abs(odd)):.1f} meV at t={int(np.argmax(abs(odd)))}")
    print(f"  max |even| (uniform)      = {np.max(abs(even)):.1f} meV")
    allpk = np.max(np.abs(S), axis=0)
    print("  peak |s| per effective band: " + " ".join(f"{x:.0f}" for x in allpk))
    if "--dump" in a:
        np.savez(a[a.index("--dump") + 1], eff=eff, S=S, odd=odd, even=even, tot=tot,
                 lo=lo, hi=hi)


if __name__ == "__main__":
    main()
