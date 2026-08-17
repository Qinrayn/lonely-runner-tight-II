# -*- coding: utf-8 -*-
"""Extend k=3 verification to r=16..25 (n up to 100)."""
import sys, math, os
_HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _HERE)
sys.path.insert(0, os.path.join(_HERE, 'lrc'))
from math import gcd
from fractions import Fraction as Fr
from ml import ml_exact, ml_float
from itertools import combinations

def check_k3_r_range(r_lo, r_hi):
    mults = [(3,3,3), (3,3,4), (3,4,3), (4,3,3), (3,4,4), (4,3,4), (4,4,3), (4,4,4)]
    counterexamples = []
    total = 0
    for r1 in range(r_lo, r_hi + 1):
        n_lo = 2 * r1 + 1
        n_hi = 4 * r1
        for n in range(n_lo, n_hi + 1):
            mid_r = [r for r in range(3, n // 2 + 1) if 2 * r <= n - 1]
            if len(mid_r) < 3 or r1 not in mid_r:
                continue
            for R in combinations(mid_r, 3):
                if r1 not in R:
                    continue
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
                    total += 1
                    val, t = ml_float(V)
                    if val >= 1.0 / n - 1e-6:
                        continue
                    ex, _ = ml_exact(V)
                    if ex < Fr(1, n):
                        counterexamples.append((n, R, W, ex))
                        print('  COUNTEREXAMPLE: n=%d R=%s W=%s LR=%s' % (n, R, W, ex))
        print('  r1=%d: done (%d total)' % (r1, total))
    return counterexamples, total

print('EXTENDED k=3 verification: r=16..25')
print('=' * 60)
ce, total = check_k3_r_range(16, 25)
print()
print('Total: %d cases, %d counterexamples' % (total, len(ce)))
if not ce:
    print('NO COUNTEREXAMPLES!')
