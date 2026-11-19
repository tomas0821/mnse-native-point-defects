#!/usr/bin/env python3
"""Relaxed version of the magnetic-order test: V_Mn and Se_Mn in the A-type and stripe AFMs-alpha orders.

Same 64-atom 4 x 2 x 2 cell, lattice, grid and cutoff as gen_stripe_defects.py, so the unrelaxed shifts can be compared
directly. Atoms relax, the lattice stays at the relaxed A-type values. The stripe host relaxes too, because it is not
force-free at the A-type geometry; the A-type host is force-free, so sd_Atype_host.out serves as its energy.
Se_Mn starts 0.60 A along -c, towards the basal Se plane, because the production relaxation found the antisite 0.77 A
that way. The A-type Se_Mn is fixed at M = -5: the unconstrained SCF creeps there slowly (47 iterations, M = -5.02), and
its energy agreed with M = -5 to 3 microeV.

Usage: gen_stripe_relax.py   (writes sr_*.in)
"""
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_stripe_defects import CELL, HERE, build  # noqa: E402

TEMPLATE = """&CONTROL
    calculation = 'relax'
    prefix      = '{name}'
    pseudo_dir  = '/home/tomas.rojas_s/pseudo'
    outdir      = './tmp_{name}'
    tprnfor     = .true.
    tstress     = .false.
    forc_conv_thr = 1.0d-4
    etot_conv_thr = 1.0d-5
    nstep       = 200
    max_seconds = 82800
/
&SYSTEM
    ibrav       = 0
    nat         = {nat}
    ntyp        = 3
    ecutwfc     = 70.0
    ecutrho     = 560.0
    occupations = 'smearing'
    smearing    = 'gaussian'
    degauss     = 0.01
    nspin       = 2{constraint}
    starting_magnetization(1) = +0.5
    starting_magnetization(2) = -0.5
    starting_magnetization(3) = +0.0
/
&ELECTRONS
    conv_thr        = 1.0d-8
    mixing_beta     = 0.1
    mixing_mode     = 'local-TF'
    electron_maxstep = 300
/
&IONS
    ion_dynamics = 'bfgs'
/
ATOMIC_SPECIES
  Mn1 54.938044 Mn.pbe-spn-kjpaw_psl.0.3.1.UPF
  Mn2 54.938044 Mn.pbe-spn-kjpaw_psl.0.3.1.UPF
  Se 78.971 Se.pbe-dn-kjpaw_psl.1.0.0.UPF
HUBBARD (ortho-atomic)
U Mn1-3d 3.7574
U Mn2-3d 3.7574
CELL_PARAMETERS angstrom
{cell}
ATOMIC_POSITIONS angstrom
{pos}
K_POINTS automatic
  2 4 3 0 0 0
"""


def write(name, atoms, constraint=""):
    lab = {+1: "Mn1", -1: "Mn2", 0: "Se"}
    pos = "\n".join(f"  {lab[a['s']]} {a['pos'][0]:.10f} {a['pos'][1]:.10f} {a['pos'][2]:.10f}" for a in atoms)
    cell = "\n".join("  " + " ".join(f"{x:.10f}" for x in v) for v in CELL)
    with open(os.path.join(HERE, name + ".in"), "w") as f:
        f.write(TEMPLATE.format(name=name, nat=len(atoms), cell=cell, pos=pos, constraint=constraint))


def main():
    vmn = ("Mn", 0, 0, 0, 0)
    for order in ("Atype", "AFMsa"):
        atoms = build(order)
        if order == "AFMsa":
            write("sr_AFMsa_host", atoms)
        write(f"sr_{order}_VMn", [a for a in atoms if a["key"] != vmn])
        sm = []
        for a in atoms:
            if a["key"] == vmn:
                b = dict(a, el="Se", s=0)
                b["pos"] = a["pos"] + [0.0, 0.0, -0.60]
                sm.append(b)
            else:
                sm.append(a)
        write(f"sr_{order}_SeMn", sm, "\n    tot_magnetization = -5" if order == "Atype" else "")


if __name__ == "__main__":
    main()
