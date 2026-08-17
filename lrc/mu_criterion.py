# -*- coding: utf-8 -*-
"""
Exact interval calculus for the Doubling Criterion (Proposition 5).

For n and m with n/2 < m < n put
    W(n,m) = { t in [0,1) : ||i t|| > 1/n  for all i in [n-1] \ {m} }.
Since the baseline [n-1] is tight, W is exactly the region where runner m alone
keeps the cover, hence ||m t|| <= 1/n on W.  Proposition 5:

    D(n,m) = ([n-1] \ {m}) u {2m}  is tight   <=>   sup_{t in W} ||m t|| <= 1/(2n).

We compute W exactly as a finite union of open intervals with rational endpoints
and evaluate mu(n,m) = max_{closure(W)} ||m t|| in exact arithmetic.  The
normalised quantity  rho(n,m) = 2 n * mu(n,m)  satisfies  rho <= 2  always, and
tightness is equivalent to rho <= 1.
"""
import os, sys, json
from fractions import Fraction as F

_HERE = os.path.dirname(os.path.abspath(__file__))


def good_set(i, n):
    """{t in [0,1) : ||i t|| > 1/n} as a list of disjoint open intervals (a,b)."""
    eps = F(1, n * i)
    out = []
    for j in range(i):
        a = F(j, i) + eps
        b = F(j + 1, i) - eps
        if a < b:
            out.append((a, b))
    return out


def intersect(A, B):
    """Intersect two sorted lists of disjoint open intervals."""
    out, ia, ib = [], 0, 0
    while ia < len(A) and ib < len(B):
        a1, a2 = A[ia]
        b1, b2 = B[ib]
        lo, hi = max(a1, b1), min(a2, b2)
        if lo < hi:
            out.append((lo, hi))
        if a2 < b2:
            ia += 1
        else:
            ib += 1
    return out


def W_intervals(n, m):
    """W(n,m) as a list of disjoint open intervals."""
    cur = [(F(0), F(1))]
    # intersect in increasing i: small i cut the circle coarsely first
    for i in range(1, n):
        if i == m:
            continue
        cur = intersect(cur, good_set(i, n))
        if not cur:
            return []
    return cur


def dist_int(x: F) -> F:
    fx = x - (x.numerator // x.denominator)
    return min(fx, 1 - fx)


def max_norm_on(a: F, b: F, m: int):
    """max of ||m t|| over the closed interval [a,b], exact."""
    best = max(dist_int(m * a), dist_int(m * b))
    # interior peaks of ||m t|| sit at t = (2j+1)/(2m)
    lo = (2 * m * a - 1) / 2      # solve (2j+1)/(2m) >= a  =>  j >= (2ma-1)/2
    import math
    j0 = math.floor(float(lo)) - 2
    j1 = math.ceil(float((2 * m * b - 1) / 2)) + 2
    for j in range(j0, j1 + 1):
        t = F(2 * j + 1, 2 * m)
        if a <= t <= b:
            best = max(best, F(1, 2))
    return best


def mu(n, m):
    """(mu, rho, witness_interval) with rho = 2 n mu; tight <=> rho <= 1."""
    W = W_intervals(n, m)
    if not W:
        return None, None, None
    best, wit = F(0), None
    for (a, b) in W:
        val = max_norm_on(a, b, m)
        if val > best:
            best, wit = val, (a, b)
    return best, 2 * n * best, wit


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 60
    smax = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    print(" n    s    m   |W|  mu(n,m)          rho=2n*mu    tight?")
    rows = []
    for n in range(5, nmax + 1):
        for s in range(2, smax + 1):
            m = n - s
            if not (2 * m > n and 1 <= m < n):
                continue
            W = W_intervals(n, m)
            if not W:
                print(f"{n:3d}  {s:3d}  {m:3d}   0   W EMPTY (!)")
                continue
            mm, rho, wit = mu(n, m)
            tight = rho <= 1
            rows.append({"n": n, "s": s, "m": m, "nW": len(W),
                         "mu": str(mm), "rho": str(rho), "rho_f": float(rho),
                         "tight": bool(tight), "witness": [str(wit[0]), str(wit[1])]})
            flag = "  <== TIGHT" if tight else ""
            print(f"{n:3d}  {s:3d}  {m:3d}  {len(W):3d}  {str(mm):16s} "
                  f"{str(rho):12s} {float(rho):.6f}{flag}")
    with open(os.path.join(_HERE, f"mu_table_{nmax}_{smax}.json"), "w") as fh:
        json.dump(rows, fh, indent=1)
