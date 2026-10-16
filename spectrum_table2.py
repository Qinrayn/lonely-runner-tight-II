# -*- coding: utf-8 -*-
"""
Table 2 support: for r=2 and 5<=n<=12, compute LR([n-1]_{2->w}) in exact
rational arithmetic over the searched window n<=w<=W0 (=8n), and report the
smallest value strictly exceeding 1/n, together with an attaining w.

Exact arithmetic (Fraction) via lrc/ml.ml_exact, which evaluates f_V at the
finite breakpoint set T(V) (peaks (2a+1)/(2v) and pairwise crossings
b/(v_i +/- v_j)) and takes the maximum.  This is the same routine that
certifies the census of Section 8.

Run:  python v3/spectrum_table2.py
"""
import sys, os
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lrc"))
from ml import ml_exact  # noqa: E402

W0FAC = 8  # W0 = 8n, the census window used throughout the paper

def Vset(n, w):
    s = set(range(1, n)) | {w}
    s.discard(2)
    return sorted(s)

print(f"{'n':>3} {'1/n':>8} {'W0':>5} {'smallest>1/n':>14} {'at w':>5} "
      f"{'<1/(n-1)?':>10} {'in [1/n,1/(n-1)) form sigma/((n-1)sigma+1)':>14}")
print("-" * 80)
for n in range(5, 13):
    target = F(1, n)
    acc = F(1, n - 1)            # accumulation point 1/(n-1)
    W0 = W0FAC * n
    best = None                   # smallest value strictly > 1/n
    best_w = None
    for w in range(n, W0 + 1):
        if w == 2:
            continue
        V = Vset(n, w)
        val, _ = ml_exact(V)
        if val > target and (best is None or val < best):
            best, best_w = val, w
    in_iso = "yes" if best < acc else "no"
    # check sigma/((n-1)sigma+1) form
    form = "---"
    if best is not None:
        num = best.numerator
        den = best.denominator
        # solve den == (n-1)*sigma + 1, sigma = (den-1)/(n-1) integer
        if (den - 1) % (n - 1) == 0:
            sigma = (den - 1) // (n - 1)
            if F(sigma, (n - 1) * sigma + 1) == best:
                form = f"sigma={sigma}"
    print(f"{n:>3} {str(target):>8} {W0:>5} {str(best):>14} {best_w:>5} "
          f"{in_iso:>10} {form:>14}")
