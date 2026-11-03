#!/usr/bin/env python3
"""Emit uideal_fold.in: primitive (ideal u) bands at every primitive k that folds onto each
point of the 3x3x2 supercell L1-Gamma-L2 path (crystal_b, 20 pts/segment, 41 points).
Supercell crystal (k1,k2,k3) -> primitive crystal (k1/3, k2/3, k3/2); folded partners add
(n1/3, n2/3, n3/2), n1,n2 in 0..2, n3 in 0..1. Ordering: K-major, partner-minor, partner 0
is the on-path point."""
import itertools, sys
src = open(sys.argv[1]).read()            # uideal_bands.in
head = src[:src.index("K_POINTS")]
head = head.replace("outdir      = './tmp_uideal'", "outdir      = './tmp_fold'")
assert "tmp_fold" in head
L1, G, L2 = (1.5, 0.0, 1.0), (0.0, 0.0, 0.0), (0.0, 1.5, 1.0)
path = []
for a, b in ((L1, G), (G, L2)):
    for i in range(20):
        path.append(tuple(a[j] + (b[j] - a[j]) * i / 20 for j in range(3)))
path.append(L2)
assert len(path) == 41
ks = []
for K in path:
    kp = (K[0] / 3, K[1] / 3, K[2] / 2)
    for n1, n2, n3 in itertools.product(range(3), range(3), range(2)):
        ks.append((kp[0] + n1 / 3, kp[1] + n2 / 3, kp[2] + n3 / 2))
assert len(ks) == 41 * 18
with open("uideal_fold.in", "w") as f:
    f.write(head + "K_POINTS crystal\n%d\n" % len(ks))
    for k in ks:
        f.write("  %.10f %.10f %.10f 1.0\n" % k)
print("wrote", len(ks), "k-points")
