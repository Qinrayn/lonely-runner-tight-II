# -*- coding: utf-8 -*-
"""
Self-falsification suite: independent attacks on the weakest joints of the Theorem A proof.

Four independent attacks on the weakest joints of the Theorem A proof and on
Theorem B's constants.  Each attack uses a code path different from the one
that produced the claim.

(A) BRUTE mid-range scan WITHOUT the effective bound: r in 3..8, all
    n in [2r+1, 6r], w in [n, 30n] -- if Theorem 1.6's bound (or its
    implementation via I) were wrong, tight sets could hide beyond it.
(B) phi-lemma attack: the proof's step 4 claims  phi(r) >= 4  ==>  2s-I >= 5
    for every mid-range n.  Verify directly for r <= 60.
(C) threshold-chain attack: step 3 claims  2s-I >= 5 ==> 4rI/(2s-I) < r(n-2).
    Verify exactly (Fractions) for r <= 60.
(D) theta-constant attack: Rosser-Schoenfeld  theta(x) > 0.84 x (x >= 101)
    and theta(x) < 1.01624 x -- verify numerically for x <= 2*10^6.
"""
import os, sys, time, math
from math import gcd
from fractions import Fraction as Fr
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LRC = os.path.join(HERE, "lrc")
sys.path.insert(0, LRC)
sys.path.insert(0, HERE)
from ml import ml_exact
from twoswap_hunt import certifies_not_tight
from allcomp_filter import Ival


def phi(r):
    return sum(1 for x in range(1, r + 1) if gcd(x, r) == 1)


def attack_A(rmax=8, wfac=30):
    print(f"(A) brute mid-range, NO effective bound: r=3..{rmax}, w <= {wfac}n")
    found = []
    t0 = time.time()
    for r in range(3, rmax + 1):
        for n in range(2 * r + 1, 6 * r + 1):
            rest = [x for x in range(1, n) if x != r]
            thr = 1.0 / n + 1e-9
            for w in range(n, wfac * n + 1):
                V = sorted(rest + [w])
                Va = np.array(V, dtype=np.int64)
                if certifies_not_tight(Va, n - 1, thr):
                    continue
                ex, _ = ml_exact(tuple(V))
                if ex == Fr(1, n):
                    found.append((r, n, w))
                    print(f"    !!! TIGHT BEYOND BOUND? r={r} n={n} w={w}")
    print(f"    tight sets found: {found if found else 'NONE'}  "
          f"({time.time()-t0:.0f}s)")
    return not found


def attack_B(rmax=60):
    print(f"(B) phi-lemma: phi(r)>=4 ==> 2s-I>=5, r<= {rmax}")
    viol = []
    for r in range(3, rmax + 1):
        if phi(r) < 4:
            continue
        for n in range(2 * r + 1, 6 * r + 1):
            s = n - r
            I = Ival(n, r)
            if I is not None and 2 * s - I < 5:
                viol.append((r, n, s, I))
    print(f"    violations: {viol if viol else 'NONE'}")
    # also: where DOES 2s-I<=4 occur?  must be exactly phi(r)<=3 territory
    occ = []
    for r in range(3, rmax + 1):
        for n in range(2 * r + 1, 6 * r + 1):
            s = n - r
            I = Ival(n, r)
            if I is not None and 2 * s - I <= 4:
                occ.append((r, n))
    rset = sorted(set(x[0] for x in occ))
    print(f"    2s-I<=4 occurs only for r in {rset} (phi: "
          f"{[phi(x) for x in rset]})")
    return not viol


def attack_C(rmax=60):
    print(f"(C) threshold chain: 2s-I>=5 ==> 4rI/(2s-I) < r(n-2), r<={rmax}")
    viol = []
    for r in range(3, rmax + 1):
        for n in range(2 * r + 1, 6 * r + 1):
            s = n - r
            I = Ival(n, r)
            if I is None or 2 * s - I < 5:
                continue
            if Fr(4 * r * I, 2 * s - I) >= r * (n - 2):
                viol.append((r, n, s, I))
    print(f"    violations: {viol if viol else 'NONE'}")
    return not viol


def attack_D(xmax=2_000_000):
    print(f"(D) theta constants, x <= {xmax}")
    sieve = np.ones(xmax + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(xmax ** 0.5) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    logs = np.zeros(xmax + 1)
    logs[sieve] = np.log(np.nonzero(sieve)[0].astype(float))
    theta = np.cumsum(logs)
    x = np.arange(101, xmax + 1, dtype=float)
    lo_ok = np.all(theta[101:] > 0.84 * x)
    hi_ok = np.all(theta[1:] < 1.01624 * np.arange(1, xmax + 1, dtype=float))
    worst_lo = (theta[101:] / x).min()
    worst_hi = (theta[1:] / np.arange(1, xmax + 1, dtype=float)).max()
    print(f"    theta(x) > 0.84x for x in [101,{xmax}]: {lo_ok} "
          f"(min ratio {worst_lo:.5f})")
    print(f"    theta(x) < 1.01624x on [1,{xmax}]: {hi_ok} "
          f"(max ratio {worst_hi:.5f})")
    return bool(lo_ok and hi_ok)


if __name__ == "__main__":
    results = {}
    results["B phi-lemma"] = attack_B()
    results["C threshold"] = attack_C()
    results["D theta"] = attack_D()
    results["A brute"] = attack_A()
    print()
    print("SELF-FALSIFICATION VERDICT")
    for k, v in results.items():
        print(f"  {'SURVIVED' if v else 'FALSIFIED'}  {k}")
