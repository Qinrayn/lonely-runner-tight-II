# -*- coding: utf-8 -*-
"""
The primorial construction -- how fast can max V grow?

Goddyn-Wong Theorem 2.3: [n-1]_{r -> mr} is tight iff gcd(r,b) > 1 for every
b in [s, ms-1], s = n-r.  In particular every PRIME p in [s, ms-1] must divide
r.  Writing theta for Chebyshev's function, this forces theta(ms-1)-theta(s-1)
<= log r, so  w = mr <= (1/2 + o(1)) n log n  (optimum at s = 2).

Conversely, taking r = primorial p# and s = 2 realises this growth: the tight
instance [r+1]_{r -> mr} with the largest admissible m has

    m = (q-1)/2,  q = smallest prime NOT dividing r,  w = mr.

This script builds the table (lower-bound data for Theorem B), verifies the GW
criterion arithmetically, and re-certifies the smaller instances by the exact
breakpoint computation and the branch-and-bound method C of paper 1.
"""
import os, sys, time
from math import gcd, log
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
LRC = os.path.join(HERE, "lrc")
sys.path.insert(0, LRC)
from ml import ml_exact
from verify_independent import prove_LR_le, f_min

PRIMES = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31]


def gw_tight(r, m, s):
    """Goddyn-Wong criterion for [r+s-1+... wait: n = r+s]: b in [s, ms-1]."""
    return all(gcd(r, b) > 1 for b in range(s, m * s))


def max_m(r, s=2):
    m = 1
    while gw_tight(r, m + 1, s):
        m += 1
    return m


def main():
    print("s = 2 primorial family: r = p#, n = r + 2, w = m r")
    print(f"{'r=p#':>12} {'n':>12} {'m':>3} {'w=mr':>14} {'w/n':>7} "
          f"{'w/(n ln n)':>11}")
    rows = []
    r = 1
    for p in PRIMES[:8]:
        r *= p
        n = r + 2
        m = max_m(r, 2)
        w = m * r
        rows.append((r, n, m, w))
        print(f"{r:>12} {n:>12} {m:>3} {w:>14} {w/n:>7.2f} "
              f"{w/(n*log(n)):>11.3f}")

    print()
    print("verification layer 1: GW criterion re-checked b by b -- implicit above")
    print("verification layer 2: exact breakpoint LR for the small members")
    for (r, n, m, w) in rows[:2]:            # r=2: n=4? r=2 -> n=4 too small; handled below
        pass
    # exact check where feasible: r=6 (n=8, w=12) and r=30 (n=32, w=90)
    for (r, m) in ((6, 2), (30, 3)):
        n = r + 2
        V = tuple(sorted([x for x in range(1, n) if x != r] + [m * r]))
        ex, t = ml_exact(V)
        print(f"  r={r:4d} n={n:4d} w={m*r:5d}: exact LR = {ex} "
              f"(=1/{n}? {ex == F(1, n)})")

    print("verification layer 3: branch-and-bound certification (method C)")
    for (r, m) in ((210, 5),):
        n = r + 2
        V = tuple(sorted([x for x in range(1, n) if x != r] + [m * r]))
        at = f_min(V, F(1, n))
        t0 = time.time()
        ok, info = prove_LR_le(V, n, max_cells=3_000_000)
        print(f"  r={r} n={n} w={m*r}: f(1/n)={at}, LR<=1/n certified={ok} "
              f"(cells={info}, {time.time()-t0:.0f}s)"
              f"  -> TIGHT" if ok and at == F(1, n) else "  -> check failed")

    print()
    print("admissible m equals (q-1)/2, q = least prime not dividing r?")
    r = 1
    for p in PRIMES[:8]:
        r *= p
        q = next(x for x in PRIMES + [37, 41, 43] if r % x != 0)
        m = max_m(r, 2)
        print(f"  r={r:>12}: m={m:>3}, (q-1)/2={(q-1)//2:>3}, equal={m == (q-1)//2}")


if __name__ == "__main__":
    main()
