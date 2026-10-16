# -*- coding: utf-8 -*-
"""
Line-by-line numeric audit of the Theorem A proof.

The proof skeleton claims, for the pair of pieces around t0 = 1/r (p = 1),
with c = the signed residue of w mod r in (-r/2, r/2], F+ = F(1), F- = F(r-1):

  (R) the right piece fits a bad interval of w
        <=>  c in [ -w/(2n) - r/n ,  r/n - w s/(n F-) ]      ... j' = 0
      and the alternative j' >= 1 is impossible when w < r(n-2);
  (L) the left piece fits
        <=>  c in [ w s/(n F+) - r/n ,  w/(2n) + r/n ]        ... j' = 0
      and j' <= -1 impossible when w < r(n-2);
  (S) signs:  r/n - w s/(n F-) < 0  and  w s/(n F+) - r/n > 0  whenever
      w >= n and s > r  -- hence (R) forces c < 0, (L) forces c > 0,
      and no w < r(n-2) survives both.

This script verifies, EXACTLY over all mid-range (n, r, w) with r in RS and
w in [n, min(4rI/(2s-I), r(n-2)-1)]:
  1. predicted (R)/(L) conditions == the exact exists-integer-j test;
  2. the sign claims (S);
  3. no w satisfies both (R) and (L)  [the contradiction].
Also: full-filter sanity for r = 2 -- the genuine tight pairs (5,7), (6,9)
MUST survive the all-components filter (necessity check of the filter itself).
"""
import os, sys, math
from math import gcd
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from allcomp_filter import Frep, Ival, pieces, filter_nr


def exact_fit(piece, w, n):
    """Exact exists-j test for one piece [a,b]."""
    a, b = piece[0], piece[1]
    lo = w * b - Fr(1, n)
    hi = w * a + Fr(1, n)
    return math.ceil(lo) <= hi


def audit(rs):
    tested = mismatches = contradictions = 0
    sign_fail = 0
    for r in rs:
        for n in range(2 * r + 1, 6 * r + 1):
            s = n - r
            I = Ival(n, r)
            if I is None or 2 * s - I <= 0:
                continue
            wmax = min((4 * r * I) // (2 * s - I), r * (n - 2) - 1)
            if wmax < n:
                continue
            F1 = Frep(n, r, 1)
            Fm = Frep(n, r, r - 1)
            pcs = [pc for pc in pieces(n, r) if pc[2] == 1]
            left = next(pc for pc in pcs if pc[3] == 'L')
            right = next(pc for pc in pcs if pc[3] == 'R')
            for w in range(n, wmax + 1):
                tested += 1
                c = w % r
                if 2 * c > r:
                    c -= r
                # predicted windows
                Rlo, Rhi = -Fr(w, 2 * n) - Fr(r, n), Fr(r, n) - Fr(w * s, n * Fm)
                Llo, Lhi = Fr(w * s, n * F1) - Fr(r, n), Fr(w, 2 * n) + Fr(r, n)
                predR = (Rlo <= c <= Rhi)
                predL = (Llo <= c <= Lhi)
                # sign claims
                if not (Rhi < 0 and Llo > 0):
                    sign_fail += 1
                # exact
                exR = exact_fit(right, w, n)
                exL = exact_fit(left, w, n)
                if predR != exR or predL != exL:
                    mismatches += 1
                    if mismatches <= 5:
                        print(f"  MISMATCH r={r} n={n} w={w} c={c}: "
                              f"pred(R,L)=({predR},{predL}) exact=({exR},{exL})")
                if exR and exL:
                    contradictions += 1
                    print(f"  BOTH COVERED (breaks proof!) r={r} n={n} w={w}")
    print(f"audit: {tested} (n,r,w) triples")
    print(f"  window-formula mismatches : {mismatches}")
    print(f"  sign-claim failures       : {sign_fail}")
    print(f"  both-pieces-covered cases : {contradictions}")
    return mismatches == 0 and sign_fail == 0 and contradictions == 0


if __name__ == "__main__":
    print("=== 1. audit of the p=1 window formulas, r=3..14 ===")
    ok = audit(range(3, 15))
    print("  => PROOF-SKELETON ARITHMETIC:", "VERIFIED" if ok else "BROKEN")
    print()
    print("=== 2. r=2 filter necessity sanity ===")
    for n in range(5, 13):
        cand, surv, ks = filter_nr(n, 2)
        print(f"  n={n:2d}: survivors = {surv if surv else 'none'}")
    print("  (5,7) and (6,9) must appear above; extras are fine -- the filter")
    print("  is necessary, not sufficient.")
