#!/usr/bin/env python3
"""Defect-cell bands inputs for Stage 5 unfolding: same k-path/nbnd as the pristine
wzsc_bands.in, geometry and &SYSTEM taken verbatim from each defect SCF input (so the
Mn_Se tot_magnetization = 3 constraint is carried over). outdir -> ./tmp_bands_<d>."""
import re
kpath = open("wzsc_bands.in").read()
kpath = kpath[kpath.index("K_POINTS"):]
for d in ("V_Mn", "V_Se", "Se_Mn", "Mn_Se"):
    s = open(f"wzsc_{d}_scf.in").read()
    s = s.replace("calculation = 'scf'", "calculation = 'bands'\n    verbosity   = 'high'")
    s = re.sub(r"\n\s*tprnfor\s*=.*", "", s)
    s = re.sub(r"\n\s*tstress\s*=.*", "", s)
    s = s.replace(f"outdir      = './tmp_wzsc_{d}'", f"outdir      = './tmp_bands_{d}'")
    s = s.replace("&SYSTEM\n", "&SYSTEM\n    nbnd        = 590\n", 1)
    s = s[:s.index("K_POINTS")] + kpath
    assert f"tmp_bands_{d}" in s and "'bands'" in s and "nbnd        = 590" in s
    c = s[s.index("&CONTROL"):s.index("/", s.index("&CONTROL"))]
    assert "verbosity" in c, "verbosity must be in &CONTROL"
    open(f"wzsc_{d}_bands.in", "w").write(s)
    print(d, "ok")
