#!/usr/bin/env python3
"""Manuscript Table 3: unfolded spin splitting of band 27 per defect cell. Primary quantities are
evaluated AT the pristine peak (|t| = 13), not as a maximum along the path (a maximum is biased high
by tracking jumps). Near Gamma (|t| <= 3) bands 27 and 28 approach each other and the tracking is unreliable, so nothing
is quoted there. Ranges: Gaussian widths 0.02/0.03/0.04/0.05 eV. Reference: the primitive cell at the
RELAXED internal parameter u = 0.3778 (bands/b_wz_bands.out, 41-point subset of the same path), because
the far field of every relaxed defect cell has left u = 0.375 and sits at u = 0.377-0.379; the unfolded
pristine supercell (u = 0.375, 239 meV) validates the method and is listed for that purpose.
Writes ../../manuscript/tab_stage5.tex; prints a plain summary."""
import os
import numpy as np
import peak_track as PT

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "..", "data", "tab_stage5.tex")
TEX = {"pristine": "Pristine", "V_Mn": r"$V_{\mathrm{Mn}}$", "V_Se": r"$V_{\mathrm{Se}}$",
       "Se_Mn": r"$\mathrm{Se}_{\mathrm{Mn}}$", "Mn_Se": r"$\mathrm{Mn}_{\mathrm{Se}}$"}
MAG = {"pristine": "0", "V_Mn": "5.25", "V_Se": "0.00", "Se_Mn": "5.00", "Mn_Se": "3$^a$"}
SIG = (0.02, 0.03, 0.04, 0.05)
T0 = 13
OK = np.abs(np.arange(41) - 20) >= 4   # path points outside the Gamma cluster
P = PT.compute(os.path.join(HERE, "w_pristine.npz"), os.path.join(HERE, "wzsc_bands.out"))
B = PT.parse(os.path.join(HERE, "..", "bands", "b_wz_bands.out"))
SUB = np.arange(0, 121, 3)
REL = (np.array(B["up"][1])[SUB, 26] - np.array(B["dn"][1])[SUB, 26]) * 1e3   # band 27, relaxed u
_t = np.arange(0, 21)
REL_ODD = (REL[20 - _t] - REL[20 + _t]) / 2
ODD_OK = _t >= 4
np.savez(os.path.join(HERE, "s_relaxed_u.npz"), S=REL, odd=REL_ODD, even=(REL[20 - _t] + REL[20 + _t]) / 2)
ref = abs(REL_ODD[T0])
rows = [f"Pristine, $u = 0.378$ & 0 & {ref:.0f} & ${REL[20 - T0]:+.0f}$/${REL[20 + T0]:+.0f}$ & 0 & 0 & -- & 1.00 \\\\"]
print(f"relaxed-u primitive: odd(t=13) {REL_ODD[T0]:.1f} meV, peak |S| {np.abs(REL).max():.1f}")
for d in ("pristine", "V_Mn", "V_Se", "Se_Mn", "Mn_Se"):
    wf = "w_pristine.npz" if d == "pristine" else f"w_{d}.npz"
    bf = "wzsc_bands.out" if d == "pristine" else f"wzsc_{d}_bands.out"
    if not (os.path.exists(os.path.join(HERE, wf)) and os.path.exists(os.path.join(HERE, bf))):
        continue
    R = {s: PT.compute(os.path.join(HERE, wf), os.path.join(HERE, bf), 27, 0.20, s) for s in SIG}
    od = [abs(float(R[s]["odd"][T0])) for s in SIG]
    ev = [abs(float(R[s]["even"][T0])) for s in SIG]
    L = [abs(float(R[s]["S"][0])) for s in SIG]
    r3 = R[0.03]
    arm = (float(r3["S"][20 - T0]), float(r3["S"][20 + T0]))
    base = P["S"] if d == "pristine" else REL
    rms = float(np.sqrt(np.nanmean((r3["S"][OK] - base[OK]) ** 2)))
    mx = float(np.nanmax(np.abs(r3["S"][OK] - base[OK])))
    rat = [abs(float(R[s]["odd"][T0])) / ref - 1 for s in SIG]
    ints = []
    for s in SIG:
        o = np.abs(R[s]["odd"])
        keep = ODD_OK & np.isfinite(o)
        ints.append(float(o[keep].sum() / np.abs(REL_ODD[keep]).sum()))
    lost = {s: int((~np.isfinite(R[s]["S"])).sum()) for s in SIG}
    np.savez(os.path.join(HERE, f"spread_{d}.npz"), odd=np.array(od), even=np.array(ev), sig=np.array(SIG))
    wm = float(np.nanmedian(r3["wgt"][:, OK]))
    wn = float(np.nanmin(r3["wgt"][:, OK]))
    fr = lambda v: (f"{min(v):.0f}" if abs(max(v) - min(v)) < 0.5 else f"{min(v):.0f}--{max(v):.0f}")
    arms = "$0$" if d == "pristine" and False else f"${arm[0]:+.0f}$/${arm[1]:+.0f}$"
    rows.append(f"{TEX[d] if d != "pristine" else "Pristine, $u = 0.375$"} & {MAG[d]} & {fr(od)} & {arms} & {fr(ev)} & {fr(L)} & {rms:.0f} & {wm:.2f} \\\\")
    print(f"{d:8s} odd(t=13) {min(od):.0f}-{max(od):.0f} | arms {arm[0]:+.0f}/{arm[1]:+.0f} | even(t=13) {min(ev):.0f}-{max(ev):.0f} "
          f"| |s(L)| {min(L):.0f}-{max(L):.0f} | rms {rms:.0f} max {mx:.0f} | w median {wm:.2f}, min {wn:.2f}"
          f" | peak vs relaxed {min(rat)*100:+.1f} to {max(rat)*100:+.1f}% | path integral {min(ints):.2f}-{max(ints):.2f} | lost k {lost}")
tex = r"""\begin{table}[t]
\centering
\caption{Unfolded spin splitting of primitive band 27 of A-type w-MnSe in the 72-atom cell (meV), at the
point of the path where the pristine splitting peaks: part odd between the arms $s_{\mathrm{odd}}$,
the signed splitting on the two arms, and the even part $s_{\mathrm{even}}$; also the
splitting at L$_{1,2}$, the rms deviation from the pristine curve at $u = 0.378$ (at $u = 0.375$ for the
pristine supercell itself), and the median spectral weight $w$ of the tracked peak (1 for a coherent
band), both over the path outside the region around $\Gamma$. Ranges: Gaussian widths 0.02--0.05~eV;
arms, rms and $w$ at 0.03~eV. At 0.02~eV the track is lost at 2 ($\mathrm{Se}_{\mathrm{Mn}}$) and 4
($\mathrm{Mn}_{\mathrm{Se}}$) of 41 points. $|M|$: net moment of the cell ($\mu_{\mathrm{B}}$). $^a$Constrained.}
\label{tab:stage5}
\small\setlength{\tabcolsep}{4pt}
\begin{tabular}{lrrrrrrr}
\toprule
Cell & $|M|$ & $s_{\mathrm{odd}}$ & Arms & $|s_{\mathrm{even}}|$ & $|s(\mathrm{L}_{1,2})|$ & rms & $w$ \\
\midrule
""" + "\n".join(rows) + r"""
\bottomrule
\end{tabular}
\end{table}
"""
open(OUT, "w").write(tex)
print("wrote", os.path.abspath(OUT))
