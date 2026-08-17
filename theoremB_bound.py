# -*- coding: utf-8 -*-
"""
Theorem B -- the extremal speed of tight single swaps.

Define  M(n) = max { max V : V a tight single-speed modification of the
baseline with n runners },  M(n) = -inf if none exists.

By Theorem A every such V with speeds > n-1 inserted is either a sporadic
((5,2,7),(6,2,9): speeds 7, 9) or has r > (n-1)/2, w = m r (m >= 2) with the
Goddyn-Wong criterion:  gcd(r,b) > 1 for all b in [s, ms-1],  s = n - r.

UPPER BOUND.  Every prime in [s, ms-1] divides r, so with P their product,
    theta(ms-1) - theta(s-1) = log P <= log r < log n .
Rosser-Schoenfeld:  theta(x) > 0.84 x  (x >= 101)  and  theta(x) < 1.01624 x.
Hence for ms-1 >= 101:
    0.84 (ms-1) - 1.01624 (s-1) < log n
    ms < 1.1905 log n + 1.2098 (s-1) + 1 + 1.2098   [absorb]
and w = m r < (ms) n / s  gives, worst at s = 2,
    w < n ( 0.5953 log n + C )         with C explicit;
for ms - 1 < 101 trivially w < 102 n / s <= 51 n.  Asymptotically (theta ~ x):
    M(n) <= (1/2 + o(1)) n log n .

LOWER BOUND.  n = p# + 2 (primorial), r = p#, s = 2, m = (q-1)/2 with q the
next prime: w = m r = ((q-1)/2) r ~ (1/2) n log n  by PNT.

This script computes M(n) exactly from the GW criterion for a large range of n
and checks it against the explicit bound  0.60 n log n + 52 n.
"""
import os, sys, time
from math import gcd, log

HERE = os.path.dirname(os.path.abspath(__file__))


def coprime_free(r, lo, hi):
    """Every b in [lo, hi] shares a factor with r?  (abort at first coprime)"""
    for b in range(lo, hi + 1):
        if gcd(b, r) == 1:
            return False
    return True


def max_m(r, s):
    """Largest m >= 1 with GW criterion; m bounded by 1 + r/(2s) (GW Lem 2.1)."""
    m = 1
    while coprime_free(r, s, (m + 1) * s - 1):
        m += 1
    return m


def M_of_n(n, s_cap=None):
    """Exact M(n) over single swaps (excluding the two sporadics, added by hand).
    s_cap limits s for speed on large n (records live at tiny s)."""
    best = None
    smax = (n - 1) // 2 if s_cap is None else min(s_cap, (n - 1) // 2)
    for s in range(2, smax + 1):          # s=1 never tight (b=1 coprime to all)
        r = n - s
        if 2 * r <= n - 1:
            break
        if not coprime_free(r, s, 2 * s - 1):   # m>=2 needs this
            continue
        m = max_m(r, s)
        if m >= 2:
            w = m * r
            if best is None or w > best[0]:
                best = (w, r, s, m)
    return best


def main():
    t0 = time.time()
    NFULL = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
    NCAP = int(sys.argv[2]) if len(sys.argv) > 2 else 40000
    bad = []
    records = []
    top_ratio = 0.0
    for n in range(5, NCAP + 1):
        cap = None if n <= NFULL else 64
        b = M_of_n(n, cap)
        if b is None:
            continue
        w, r, s, m = b
        bound = 0.60 * n * log(n) + 52 * n
        if w > bound:
            bad.append((n, b))
        ratio = w / (n * log(n))
        if ratio > top_ratio + 1e-12:
            top_ratio = ratio
            records.append((n, w, r, s, m, round(ratio, 4)))
    print(f"scanned n <= {NCAP} (full s for n <= {NFULL}, s <= 64 beyond); "
          f"{time.time()-t0:.0f}s")
    print(f"violations of  w <= 0.60 n ln n + 52 n : "
          f"{bad if bad else 'NONE'}")
    print("record ratios w / (n ln n)  [n, w, r, s, m, ratio]:")
    for rec in records[-12:]:
        print("   ", rec)


if __name__ == "__main__":
    main()
