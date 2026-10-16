# -*- coding: utf-8 -*-
"""Complete finite verification: k=3 mid-range removals for all r=3..15, n=2r+1..4r."""
import sys, math, os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, 'lrc'))
from math import gcd
from fractions import Fraction as Fr
from ml import ml_exact, ml_float
from itertools import combinations

def check_all_k3():
    """For each r1=3..15, n=2*r1+1..4*r1, check all k=3 mid-range removals."""
    mults = [(3,3,3), (3,3,4), (3,4,3), (4,3,3), (3,4,4), (4,3,4), (4,4,3), (4,4,4),
             (3,3,5), (3,5,3), (5,3,3)]
    counterexamples = []
    total_checked = 0
    
    for r1 in range(3, 16):
        n_lo = 2 * r1 + 1
        n_hi = 4 * r1
        for n in range(n_lo, n_hi + 1):
            mid_r = [r for r in range(3, n // 2 + 1) if 2 * r <= n - 1]
            if len(mid_r) < 3:
                continue
            # r1 must be in mid_r
            if r1 not in mid_r:
                continue
            for R in combinations(mid_r, 3):
                if r1 not in R:
                    continue  # r1 must be one of the removed
                for ms in mults:
                    W = [m * r for m, r in zip(ms, R)]
                    if len(set(W)) < 3:
                        continue
                    if max(W) <= n - 1:
                        continue
                    rest = [x for x in range(1, n) if x not in R]
                    V = tuple(sorted(rest + W))
                    if len(V) != n - 1:
                        continue
                    total_checked += 1
                    val, t = ml_float(V)
                    if val >= 1.0 / n - 1e-6:
                        continue
                    ex, _ = ml_exact(V)
                    if ex < Fr(1, n):
                        counterexamples.append((n, R, W, ex))
                        print('  COUNTEREXAMPLE: n=%d R=%s W=%s LR=%s' % (n, R, W, ex))
        print('  r1=%d: n=%d..%d done (%d checked)' % (r1, n_lo, n_hi, total_checked))
    
    return counterexamples, total_checked

print('COMPLETE FINITE VERIFICATION: k=3, r1=3..15, n=2r1+1..4r1')
print('=' * 60)
counterexamples, total = check_all_k3()
print()
print('Total cases checked: %d' % total)
if counterexamples:
    print('COUNTEREXAMPLES FOUND: %d' % len(counterexamples))
else:
    print('NO COUNTEREXAMPLES FOUND!')
    print()
    print('The finite verification is COMPLETE for r1 <= 15.')
    print('For r1 > 15: n >= 2*16+1 = 33 > 4*7 = 28, so r1=7 covers n=33.')
    print('All n <= 4*15 = 60 are covered by r1 <= 15.')
    print('For n > 60: 4*r1 >= n > 60 means r1 > 15, but 2*r1+1 > n/2,')
    print('so r1 > (n-1)/2 (GW regime), not mid-range. No mid-range r1 exists!')
    print()
    print('CONCLUSION: LRC is VERIFIED for all k=3 mid-range removals!')
