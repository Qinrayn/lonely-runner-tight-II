# -*- coding: utf-8 -*-
"""Extend k>=3 mid-range verification to n=26..35."""
import sys, math, os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, 'lrc'))
from math import gcd
from fractions import Fraction as Fr
from ml import ml_exact, ml_float
from itertools import combinations

def check_k3(n, wmax_m=6):
    """Check all k=3 mid-range removals for counterexamples."""
    mid_r = [r for r in range(3, n//2+1) if 2*r <= n-1]
    if len(mid_r) < 3:
        return 'skip'
    # Multipliers to try (m>=3 needed for v > n-1 with r <= (n-1)/2)
    mults = [(3,3,3), (3,3,4), (3,4,3), (4,3,3), (3,4,4), (4,3,4), (4,4,3), (4,4,4),
             (3,3,5), (3,5,3), (5,3,3), (3,5,5), (5,3,5), (5,5,3), (5,5,5),
             (4,4,5), (4,5,4), (5,4,4), (4,5,5), (5,4,5), (5,5,4)]
    for R in combinations(mid_r, 3):
        for ms in mults:
            W = [m*r for m, r in zip(ms, R)]
            if len(set(W)) < 3:
                continue
            # At least one w > n-1
            if max(W) <= n-1:
                continue
            rest = [x for x in range(1, n) if x not in R]
            V = tuple(sorted(rest + W))
            if len(V) != n-1:
                continue
            val, t = ml_float(V)
            if val >= 1.0/n - 1e-6:
                continue
            ex, _ = ml_exact(V)
            if ex < Fr(1, n):
                return 'COUNTEREXAMPLE: R=%s W=%s LR=%s' % (R, W, ex)
    return 'OK'

print('k=3 mid-range verification: n=26..35')
for n in range(26, 36):
    result = check_k3(n)
    print('  n=%d: %s' % (n, result))
print('Done.')
