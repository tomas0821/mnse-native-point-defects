#!/usr/bin/env python3
"""Shift of the formation energies between the stripe AFMs-alpha and A-type orders of w-MnSe, unrelaxed (sd_*) and relaxed (sr_*_L*),
in the same 64-atom 4 x 2 x 2 cell, at fixed dmu_Se, each order on its own host. Relaxed A-type host = sd_Atype_host (force-free).
usage: order_shift_relaxed.py"""
import glob
import re

RY = 13.605693123


def energy(f):
    return float(re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)", open(f, errors="ignore").read(), re.M)[-1])


def last_leg(tag):
    return sorted(glob.glob(f"{tag}_L*.out"))[-1]


def ef(e_def, e_host, n_mn=-1):
    return (e_def - e_host - n_mn * e_host / 32) * RY


h = {"un": {o: energy(f"sd_{o}_host.out") for o in ("Atype", "AFMsa")},
     "re": {"Atype": energy("sd_Atype_host.out"), "AFMsa": energy(last_leg("sr_AFMsa_host"))}}
print(f"host, stripe minus A-type: unrelaxed {(h['un']['AFMsa'] - h['un']['Atype']) * RY * 1e3 / 32:+.2f}, "
      f"relaxed {(h['re']['AFMsa'] - h['re']['Atype']) * RY * 1e3 / 32:+.2f} meV/f.u.")
for d in ("VMn", "SeMn"):
    un = ef(energy(f"sd_AFMsa_{d}.out"), h["un"]["AFMsa"]) - ef(energy(f"sd_Atype_{d}.out"), h["un"]["Atype"])
    re_ = ef(energy(last_leg(f"sr_AFMsa_{d}")), h["re"]["AFMsa"]) - ef(energy(last_leg(f"sr_Atype_{d}")), h["re"]["Atype"])
    print(f"{d}: stripe minus A-type {un:+.3f} eV unrelaxed, {re_:+.3f} eV relaxed")
