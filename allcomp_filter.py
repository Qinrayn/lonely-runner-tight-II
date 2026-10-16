# -*- coding: utf-8 -*-
"""
The all-components filter (the engine for Theorem A).

Tightness of [n-1]_{r->w} demands that w alone covers EVERY component of
U(n,r), not just the longest one (paper 1 used only the longest).  In the
mid-range 2r <= n-1 the components near t0 = p/r (gcd(p,r)=1) are the two
pieces
    left :  [p/r - s/(r n F(u)),   p/r - 1/(2rn)]      u = p^{-1} mod r
    right:  [p/r + 1/(2rn),        p/r + s/(r n F(r-u))]
(nonempty iff F < 2s).  A piece [a,b] fits inside a bad interval
[j/w - 1/(wn), j/w + 1/(wn)] iff there is an integer j with
    w*b - 1/n  <=  j  <=  w*a + 1/n .
The filter: w survives iff every nonempty piece admits such a j.

Conjecture behind Theorem A: for r >= 3, no w in [n, 4rI/(2s-I)] survives.
This script tests the filter EXACTLY (rationals), reports survivors, and for
each killed w records a witness piece -- data to guide the eventual proof.
"""
import os, sys, json
from math import gcd
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))


def Frep(n, r, x):
    reps = [i for i in range(1, n) if i != r and (i - x) % r == 0]
    return max(reps) if reps else None


def pieces(n, r):
    """All nonempty component pieces of U(n,r) in the mid-range, exact."""
    s = n - r
    out = []
    for p in range(1, r):
        if gcd(p, r) != 1:
            continue
        u = pow(p, -1, r)
        for sign, x in ((-1, u), (+1, r - u)):
            Fx = Frep(n, r, x % r if x % r else r)  # x in {u, r-u}, never 0
            if Fx is None or Fx >= 2 * s:
                continue                            # piece empty
            t0 = Fr(p, r)
            inner = Fr(1, 2 * r * n)
            outer = Fr(s, r * n * Fx)
            if sign < 0:
                out.append((t0 - outer, t0 - inner, p, 'L'))
            else:
                out.append((t0 + inner, t0 + outer, p, 'R'))
    return out


def Ival(n, r):
    best = None
    for u in range(1, r):
        if gcd(u, r) != 1:
            continue
        for x in (u, r - u):
            v = Frep(n, r, x)
            if v is not None and (best is None or v < best):
                best = v
    return best


def filter_nr(n, r):
    """Return (candidate count, survivors list, kill statistics)."""
    s = n - r
    I = Ival(n, r)
    if I is None or 2 * s - I <= 0:
        return 0, [], {}
    wmax = (4 * r * I) // (2 * s - I)
    pcs = pieces(n, r)
    survivors, killstat = [], {}
    for w in range(n, wmax + 1):
        ok = True
        for (a, b, p, side) in pcs:
            lo = w * b - Fr(1, n)          # j >= lo
            hi = w * a + Fr(1, n)          # j <= hi
            # integer in [lo, hi]?
            import math
            jmin = math.ceil(lo)
            if jmin > hi:
                ok = False
                killstat[(p, side)] = killstat.get((p, side), 0) + 1
                break
        if ok:
            survivors.append(w)
    return wmax - n + 1, survivors, killstat


if __name__ == "__main__":
    rs = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 \
        else [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    total_surv = {}
    for r in rs:
        surv_r = []
        for n in range(2 * r + 1, 6 * r + 1):
            cand, surv, ks = filter_nr(n, r)
            if surv:
                surv_r.append((n, surv))
        total_surv[r] = surv_r
        flat = [w for (_, ws) in surv_r for w in ws]
        print(f"r={r:3d}: candidates over all n exist; "
              f"FILTER SURVIVORS: {surv_r if surv_r else 'NONE'}")
    print()
    any_surv = any(v for v in total_surv.values())
    print("=> all-components filter alone kills every candidate:" ,
          "NO" if any_surv else "YES  (Theorem A reduces to formalising the filter)")
