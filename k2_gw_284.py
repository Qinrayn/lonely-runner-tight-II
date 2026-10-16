# -*- coding: utf-8 -*-
"""Branch-and-bound certification of the second k-GW member (n=284).

V = [283] \\ {280, 282} U {560, 564};  claims LR(V) = 1/284.
Method C of lrc/verify_independent.py: adaptive bisection with the
Lipschitz bound, exact rational arithmetic throughout.  Also prints the
cell count, which the paper quotes.
"""
import sys
from fractions import Fraction as F

sys.path.insert(0, "lrc")
from verify_independent import prove_LR_le, f_min

n = 284
R = (280, 282)
W = (560, 564)
V = tuple(sorted([x for x in range(1, n) if x not in R] + list(W)))
assert len(V) == n - 1

at = f_min(V, F(1, n))          # f_V(1/n): speed 1 gives 1/n
ok, cells = prove_LR_le(V, n, max_cells=40_000_000)
print(f"n=284 k-GW member V=[283]\\{{280,282}}+{{560,564}}:")
print(f"  f_V(1/n) = {at}  (== 1/284: {at == F(1, n)})")
print(f"  LR(V) <= 1/284 certified: {ok}   (branch-and-bound cells: {cells})")
print("  => TIGHT" if ok and at == F(1, n) else "  => CHECK FAILED")
sys.exit(0 if (ok and at == F(1, n)) else 1)
