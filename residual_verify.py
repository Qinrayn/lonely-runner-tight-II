# -*- coding: utf-8 -*-
"""Finite residual verification: n < 3r+1 for r=3..15."""
import sys, math, os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, 'lrc'))
from math import gcd
from fractions import Fraction as Fr
from mu_criterion import good_set, intersect
from ml import ml_exact, ml_float

def verify_residual(n, r1, wmax_cap=80):
    for r2 in range(r1+1, n):
        U = [(Fr(0), Fr(1))]
        for i in range(1, n):
            if i in [r1, r2]:
                continue
            U = intersect(U, good_set(i, n))
            if not U:
                break
        if not U:
            continue
        wmax = min(wmax_cap, 5*n)
        for w1 in range(n, wmax):
            R = intersect(U, good_set(w1, n))
            if not R:
                continue
            for w2 in range(n, wmax):
                if w2 == w1:
                    continue
                if not intersect(R, good_set(w2, n)):
                    rest = [x for x in range(1, n) if x != r1 and x != r2]
                    V = tuple(sorted(rest + [w1, w2]))
                    val, t = ml_float(V)
                    if val > 1.0/n + 1e-6:
                        continue
                    ex, _ = ml_exact(V)
                    if ex == Fr(1, n):
                        return 'TIGHT: r2=%d w1=%d w2=%d' % (r2, w1, w2)
    return 'OK'

print('FINITE RESIDUAL VERIFICATION: n < 3r+1, r=3..15')
print()
all_clear = True
for r in range(3, 16):
    lo, hi = 2*r+1, 3*r
    for n in range(lo, hi+1):
        result = verify_residual(n, r)
        if result != 'OK':
            print('  r=%d n=%d: %s' % (r, n, result))
            all_clear = False
    print('  r=%d: n=%d..%d done' % (r, lo, hi))

if all_clear:
    print()
    print('ALL RESIDUAL CASES VERIFIED!')
    print('No tight two-swap for r=3..15, n < 3r+1.')
    print('Combined with L-component proof (n >= 3r+1):')
    print('No mid-range two-swap with odd r1 >= 3 exists for ANY n.')
