# -*- coding: utf-8 -*-
"""
General two-swap hunt: attack the open converse of Goddyn--Wong Theorem 3.1.

We remove a pair {r1,r2} from the baseline [n-1] and insert two arbitrary new
speeds w1 < w2 (both > n-1, bounded by WFAC*n), test tightness exactly, and
classify the outcome:

  * "aligned"  : r1|w1 and r2|w2 (or r1|w2 and r2|w1 with matching multipliers)
                 and each accelerated runner satisfies GW's per-runner GCD
                 condition  ->  predicted tight by GW Thm 3.1;
  * "CROSSED"  : the divisibility pattern is crossed (r1|w2, r2|w1 only), or
  * "VIOLATION": tight although some individual runner fails its GCD condition
                 ->  resolves GW's open question in the negative.

Generalised Proposition R is also tested: for tight instances with
2*r_i > n-1 we must have r_i | w1 or r_i | w2 for each i.
"""
import os, sys, time
from math import gcd
from itertools import combinations
from fractions import Fraction
import numpy as np
from numba import njit

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ml import ml_exact


@njit(cache=True)
def _f_at(v, k, t):
    best = 1.0
    for i in range(k):
        x = v[i] * t
        d = abs(x - round(x))
        if d < best:
            best = d
    return best


@njit(cache=True)
def certifies_not_tight(v, k, thr):
    for i in range(k):
        d = 2 * v[i]
        for num in range(1, d, 2):
            if _f_at(v, k, num / d) > thr:
                return True
    for i in range(k):
        for j in range(i + 1, k):
            for d in (v[i] + v[j], v[j] - v[i]):
                if d <= 0:
                    continue
                for num in range(1, d):
                    if _f_at(v, k, num / d) > thr:
                        return True
    return False


def gw_cond(n, r, m):
    """GW per-runner GCD condition for r -> m*r."""
    s = n - r
    return all(gcd(r, b) > 1 for b in range(s, m * s))


def analyse(n, r1, r2, w1, w2):
    """Return a verdict string for a tight instance."""
    notes = []
    for r in (r1, r2):
        divs = [w for w in (w1, w2) if w % r == 0]
        if not divs:
            notes.append(f"PropR-general FAILS for r={r}")
    # try to match each r with a w it divides
    pairings = []
    for (a, b) in ((w1, w2), (w2, w1)):
        if a % r1 == 0 and b % r2 == 0:
            pairings.append((a // r1, b // r2, a, b))
    if not pairings:
        return "CROSSED/UNPAIRED " + ";".join(notes), None
    verdicts = []
    for (m1, m2, a, b) in pairings:
        c1, c2 = gw_cond(n, r1, m1), gw_cond(n, r2, m2)
        verdicts.append((c1 and c2, m1, m2, c1, c2))
    if any(v[0] for v in verdicts):
        v = next(v for v in verdicts if v[0])
        return f"aligned+GWcond (m1={v[1]},m2={v[2]})", True
    v = verdicts[0]
    return (f"!!! GW-CONVERSE VIOLATION: m1={v[1]}(cond={v[3]}) "
            f"m2={v[2]}(cond={v[4]})"), False


def hunt(nmax, wfac=4, nmin=6):
    t0 = time.time()
    results = []
    for n in range(nmin, nmax + 1):
        thr = 1.0 / n + 1e-9
        base = list(range(1, n))
        WMAX = wfac * n
        for r1, r2 in combinations(range(2, n), 2):
            rest = [x for x in base if x not in (r1, r2)]
            for w1, w2 in combinations(range(n, WMAX + 1), 2):
                V = rest + [w1, w2]
                if len(set(V)) != n - 1:
                    continue
                Vs = np.array(sorted(V), dtype=np.int64)
                if certifies_not_tight(Vs, n - 1, thr):
                    continue
                Vt = tuple(int(x) for x in Vs)
                ex, tstar = ml_exact(Vt)
                if ex < Fraction(1, n):
                    print(f"!!!!!!!!!! LRC COUNTEREXAMPLE n={n}: {Vt} ML={ex}")
                    continue
                if ex != Fraction(1, n):
                    continue
                verdict, ok = analyse(n, r1, r2, w1, w2)
                results.append((n, r1, r2, w1, w2, verdict))
                print(f"TIGHT n={n}: remove {{{r1},{r2}}} add {{{w1},{w2}}}  -> {verdict}")
        print(f"  ... n={n} done ({time.time()-t0:.0f}s)")
    print("\n=== summary ===")
    viol = [r for r in results if "VIOLATION" in r[5] or "CROSSED" in r[5]]
    print(f"tight two-swaps found: {len(results)};  anomalous: {len(viol)}")
    for v in viol:
        print("  ", v)
    return results


if __name__ == "__main__":
    hunt(int(sys.argv[1]), int(sys.argv[2]) if len(sys.argv) > 2 else 4,
         int(sys.argv[3]) if len(sys.argv) > 3 else 6)
