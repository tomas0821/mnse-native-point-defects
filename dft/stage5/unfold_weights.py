#!/usr/bin/env python3
"""Spin-resolved Popescu-Zunger unfolding weights from QE (7.x) wfc{up,dw}<ik>.dat files.

For a supercell whose lattice is diag(M) times the primitive one (here M = 3,3,2), each
supercell K on the band path was chosen EQUAL to the primitive on-path k. The primitive
Bloch component at k = K therefore consists of the supercell plane waves K+G whose
supercell Miller indices satisfy  m_i % M_i == 0  for all i. The unfolding weight of band n
is
        w_n(K) = sum_{G in that sublattice} |c_n(G)|^2  /  sum_G |c_n(G)|^2
(pseudo-wavefunction coefficients, normalised per band -- the standard approximation for
PAW/USPP, as in BandUP). Pristine cell: w is 0 or 1. Defect cell: fractional.

Record layout written by QE io_base.f90 :: write_wfc (Fortran sequential, 4-byte markers):
  1: ik(i4) xk(3 f8) ispin(i4) gamma_only(l4) scalef(f8)
  2: ngw igwx npol nbnd (4 i4)
  3: b1 b2 b3 (9 f8)
  4: mill(3, igwx) (i4)
  5..: evc(npol*igwx) complex(c16), one record per band

usage: unfold_weights.py SAVE_DIR OUT.npz [--M 3 3 2] [--nproc N]
Output npz: w[spin, ik, band] (spin 0 = up, 1 = dn), xk[ik, 3], nbnd, M.
"""
import os
import sys
import struct
import numpy as np


def _rec(f):
    """Read one Fortran sequential record; return raw bytes."""
    h = f.read(4)
    if len(h) < 4:
        raise EOFError
    n = struct.unpack("<i", h)[0]
    data = f.read(n)
    t = struct.unpack("<i", f.read(4))[0]
    assert t == n, f"record marker mismatch {n} != {t}"
    return data


def read_weights(path, M):
    with open(path, "rb") as f:
        r1 = _rec(f)
        assert len(r1) == 44, f"unexpected header record length {len(r1)}"
        ik, = struct.unpack_from("<i", r1, 0)
        xk = np.frombuffer(r1, "<f8", 3, 4)
        ispin, gam = struct.unpack_from("<ii", r1, 28)
        assert gam == 0, "gamma_only files not supported"
        ngw, igwx, npol, nbnd = struct.unpack("<4i", _rec(f))
        assert npol == 1, "noncollinear not supported"
        _rec(f)  # b1 b2 b3
        mill = np.frombuffer(_rec(f), "<i4").reshape(igwx, 3)
        mask = np.ones(igwx, bool)
        for i in range(3):
            mask &= (mill[:, i] % M[i]) == 0
        w = np.empty(nbnd)
        for n in range(nbnd):
            c = np.frombuffer(_rec(f), "<c16")
            a = np.abs(c) ** 2
            w[n] = a[mask].sum() / a.sum()
    return ik, ispin, xk.copy(), w, int(mask.sum()), igwx


def _job(args):
    path, M = args
    return path, read_weights(path, M)


def main():
    a = sys.argv[1:]
    save, out = a[0], a[1]
    M = (3, 3, 2)
    if "--M" in a:
        i = a.index("--M")
        M = tuple(int(x) for x in a[i + 1:i + 4])
    nproc = 1
    if "--nproc" in a:
        nproc = int(a[a.index("--nproc") + 1])
    files = sorted(f for f in os.listdir(save) if f.startswith("wfc") and f.endswith(".dat"))
    nk = max(int(f[5:-4]) for f in files)
    paths = [os.path.join(save, f"wfc{s}{k}.dat") for s in ("up", "dw") for k in range(1, nk + 1)]
    missing = [p for p in paths if not os.path.exists(p)]
    assert not missing, f"missing {missing[:3]}"
    jobs = [(p, M) for p in paths]
    if nproc > 1:
        from multiprocessing import Pool
        with Pool(nproc) as pool:
            res = dict(pool.imap_unordered(_job, jobs))
    else:
        res = dict(map(_job, jobs))
    nbnd = len(next(iter(res.values()))[3])
    w = np.zeros((2, nk, nbnd))
    xk = np.zeros((nk, 3))
    for p, (ik, ispin, x, wk, nsub, igwx) in res.items():
        s = 0 if os.path.basename(p).startswith("wfcup") else 1
        w[s, ik - 1] = wk
        xk[ik - 1] = x
    np.savez(out, w=w, xk=xk, nbnd=nbnd, M=np.array(M))
    print(f"wrote {out}: nk={nk} nbnd={nbnd} M={M} "
          f"sum_w(up,k1)={w[0,0].sum():.3f} sum_w(dn,k1)={w[1,0].sum():.3f}")


if __name__ == "__main__":
    main()
