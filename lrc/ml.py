# -*- coding: utf-8 -*-
"""
Exact computation of the maximum loneliness ML(v) for the Lonely Runner problem.

    ML(v) = max_{t in [0,1)} min_i ||v_i t||,   v = distinct positive integers.

Rigor backbone: the Breakpoint Lemma -- the maximum is attained at a point of
    T(v) = { (2a+1)/(2 v_i) } U { b/(v_i+v_j) } U { b/(v_j-v_i) }  (i<j),
so evaluating f at these finitely many rationals with exact arithmetic yields
a *proof-grade* value of ML(v).
"""
from fractions import Fraction
from math import gcd
from functools import reduce
import numpy as np


def dist_to_int_frac(x: Fraction) -> Fraction:
    """||x|| for a Fraction x, exact."""
    fx = x - (x.numerator // x.denominator)  # fractional part in [0,1)
    return min(fx, 1 - fx)


def f_exact(v, t: Fraction) -> Fraction:
    """f_v(t) = min_i ||v_i t||, exact."""
    return min(dist_to_int_frac(vi * t) for vi in v)


def candidate_denominators(v):
    """All denominators from the Breakpoint Lemma (with multiplicity removed)."""
    dens = set()
    for vi in v:
        dens.add(2 * vi)          # peaks (odd numerators only, filtered later)
    k = len(v)
    for i in range(k):
        for j in range(i + 1, k):
            dens.add(v[i] + v[j])
            d = abs(v[j] - v[i])
            if d:
                dens.add(d)
    return dens


def ml_exact(v):
    """Exact ML(v) plus one maximizing t. Returns (Fraction, Fraction)."""
    v = sorted(v)
    best, best_t = Fraction(0), Fraction(0)
    # peaks: t = (2a+1)/(2 vi)
    for vi in v:
        den = 2 * vi
        for num in range(1, den, 2):
            t = Fraction(num, den)
            val = f_exact(v, t)
            if val > best:
                best, best_t = val, t
    # crossings: t = b/(vi+vj) and b/(vj-vi)
    seen = set()
    k = len(v)
    for i in range(k):
        for j in range(i + 1, k):
            for den in (v[i] + v[j], v[j] - v[i]):
                if den <= 0 or den in seen:
                    continue
                seen.add(den)
                for num in range(1, den):
                    t = Fraction(num, den)
                    val = f_exact(v, t)
                    if val > best:
                        best, best_t = val, t
    return best, best_t


def ml_float(v):
    """Fast float ML(v) via numpy over the same candidate set (prescreening)."""
    v = np.asarray(sorted(v), dtype=np.int64)
    ts = []
    for vi in v:
        den = 2 * vi
        ts.append(np.arange(1, den, 2, dtype=np.float64) / den)
    k = len(v)
    seen = set()
    for i in range(k):
        for j in range(i + 1, k):
            for den in (int(v[i] + v[j]), int(v[j] - v[i])):
                if den <= 0 or den in seen:
                    continue
                seen.add(den)
                ts.append(np.arange(1, den, dtype=np.float64) / den)
    t = np.concatenate(ts)
    # f(t) = min_i || v_i t ||
    x = np.outer(t, v)
    x = np.abs(x - np.rint(x))
    fvals = x.min(axis=1)
    idx = int(fvals.argmax())
    return float(fvals[idx]), float(t[idx])


def normalize(v):
    """Sort, take |.|, divide by gcd. ML is invariant under these."""
    v = sorted(abs(int(x)) for x in v)
    g = reduce(gcd, v)
    return tuple(x // g for x in v)


if __name__ == "__main__":
    # Sanity checks against known theorems (LRC tight baseline {1..k} -> 1/(k+1))
    for k in range(2, 7):
        v = list(range(1, k + 1))
        val, t = ml_exact(v)
        print(f"k={k}  v={v}  ML={val}  (expect 1/{k+1})  at t={t}")
    # A few random instances must satisfy ML >= 1/(k+1) since LRC holds for k<=12
    import random
    random.seed(1)
    for _ in range(5):
        k = 4
        v = normalize(random.sample(range(1, 60), k))
        val, t = ml_exact(v)
        ok = val >= Fraction(1, k + 1)
        print(f"v={v}  ML={val} ~ {float(val):.5f}  >=1/{k+1}: {ok}")
