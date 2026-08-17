# -*- coding: utf-8 -*-
"""
Theoretical exclusion of tight two-swaps R = {q, r'} containing an arbitrary
small removal q (Section 6.2 of the paper).  No restriction on r' (any
regime), no upper bound on the inserted speeds.

Chain of results (all verified here in exact rational arithmetic):

(L1') Windows at 1/q (Lemma "windows"): with F_{+-1} the largest retained
      speed in the class +-1 mod q and m0*q the smallest retained multiple
      of q (m0 in {2,3} once 3q <= n-1),
        W_- = (1/q - (1/q - 1/n)/F_1 ,  1/q - 1/(m0*q*n)),
        W_+ = (1/q + 1/(m0*q*n)      ,  1/q + (1/q - 1/n)/F_{-1})
      are contained in U(n, {q, r'}); if V is tight with speeds >= n, the
      bad sets of w1, w2 must cover W_- u W_+.
      [Verified here: they are exactly the two components of U(n, {q, r'})
      adjacent to 1/q, over the full grid n <= 60, 4q <= n-1, all r'.]

(L2') Edge reach at 1/q (Lemma "edge reach"): a bad interval of a speed
      w >= n covering points of W_+ arbitrarily close to the inner edge
      a = 1/q + 1/(m0*q*n) forces w >= m0*(n-q) >= 2(n-q); the alternative
      (q | w and w <= m0*q < n) is void.
      [Verified by brute force over all bad intervals, integer arithmetic.]

(T)  Density inequality (Theorem "small removal"): with
        ell* = (n - 2q + 1) / (2 q n (n-1)),
        ell* * (1 - 4/n) > (1 + 2/n)/((n-q) n) + 2 (1 + 2/n)/n^2        (*)
      exclusion follows.  (*) holds for every q >= 2 with n >= max(40, 10q)
      (crude lemma, verified exactly), with refined exact thresholds n0(q)
      (verified with no later failure), and the asymptotic threshold is
      gamma = (9 - sqrt(57))/12 ~ 0.1208 (cross-multiplied difference is a
      quartic in n with leading form n^4 (6 c^2 - 9 c + 1), c = q/n).

Positive control: q = 2 with parity-aware F (F = n-2 for odd n) reproduces
the thresholds of r2_density.py exactly (odd n >= 19, even n >= 20).

Output: every claim is checked; the script exits nonzero on failure.
"""
import sys
from fractions import Fraction as Fr

sys.path.insert(0, "lrc")
from mu_criterion import good_set, intersect


def U_components(n, removed):
    cur = [(Fr(0), Fr(1))]
    for i in range(1, n):
        if i in removed:
            continue
        cur = intersect(cur, good_set(i, n))
        if not cur:
            return []
    return cur


def predicted_windows(n, q, rp):
    retained = [i for i in range(1, n) if i not in (q, rp)]
    F1 = max(i for i in retained if i % q == 1)
    Fm1 = max(i for i in retained if i % q == q - 1)
    m0 = min(i for i in retained if i % q == 0) // q
    betaL = (Fr(1, q) - Fr(1, n)) / F1
    betaR = (Fr(1, q) - Fr(1, n)) / Fm1
    inner = Fr(1, m0 * q * n)
    return (Fr(1, q) - betaL, Fr(1, q) - inner), (Fr(1, q) + inner, Fr(1, q) + betaR), m0


def main():
    ok = True

    print("=" * 64)
    print("(L1') windows at 1/q are exactly the adjacent components")
    print("=" * 64)
    mismatch = checked = 0
    for n in range(9, 61):
        for q in range(2, (n - 1) // 4 + 1):
            for rp in range(1, n):
                if rp == q:
                    continue
                predL, predR, m0 = predicted_windows(n, q, rp)
                assert predL[0] < predL[1] and predR[0] < predR[1], (n, q, rp)
                comps = U_components(n, {q, rp})
                checked += 1
                if predL not in comps or predR not in comps:
                    mismatch += 1
                    if mismatch <= 5:
                        print("  MISMATCH", n, q, rp)
    ok &= (mismatch == 0)
    print(f"  instances checked: {checked}, mismatches: {mismatch}")

    print("=" * 64)
    print("(L2') edge reach: covering the inner edge forces w >= m0(n-q)")
    print("=" * 64)
    viol = tot = 0
    for n in range(9, 51):
        for q in range(2, (n - 1) // 4 + 1):
            for m0 in (2, 3):
                # interval [j/w - 1/(wn), j/w + 1/(wn)] covers points of W_+
                # arbitrarily close to a = 1/q + 1/(m0 q n)  iff  (closed)
                #   j/w - 1/(wn) <= a < j/w + 1/(wn)
                # integer form (multiply by w*q*m0*n):
                #   j*q*m0*n - q*m0 <= w*m0*n + w*q   and
                #   j*q*m0*n + q*m0 >  w*m0*n + w*q
                for w in range(n, 5 * n):
                    for j in range(0, w + 1):
                        # a = (m0*n+1)/(m0*q*n); j/w - 1/(wn) <= a < j/w + 1/(wn)
                        # multiplied by w*m0*q*n:  j*q*m0*n - q*m0 <= w*(m0*n+1)
                        # and  j*q*m0*n + q*m0 > w*(m0*n+1)
                        lhs = j * q * m0 * n
                        if lhs - q * m0 <= w * (m0 * n + 1) < lhs + q * m0:
                            tot += 1
                            if not (w >= m0 * (n - q)
                                    or (w % q == 0 and w <= m0 * q)):
                                viol += 1
                                if viol <= 5:
                                    print("  VIOLATION", n, q, m0, w, j)
    ok &= (viol == 0)
    print(f"  covering intervals: {tot}, violations: {viol}")

    print("=" * 64)
    print("(T) density inequality (*): crude lemma, exact thresholds, limit")
    print("=" * 64)

    def holds(n, q):
        ell = Fr(n - 2 * q + 1, 2 * q * n * (n - 1))
        lhs = ell * (1 - Fr(4, n))
        rhs = (1 + Fr(2, n)) / ((n - q) * n) + 2 * (1 + Fr(2, n)) / n ** 2
        return lhs > rhs

    # crude lemma: q >= 2, n >= max(40, 10q); complete for q <= 20, n <= 2000,
    # boundary band + spot-check beyond
    crude_fail = [(n, q) for q in range(2, 21) for n in range(max(40, 10 * q), 2001)
                  if not holds(n, q)]
    crude_fail += [(n, q) for q in range(21, 201) for n in range(max(40, 10 * q), max(40, 10 * q) + 40)
                   if not holds(n, q)]
    crude_fail += [(n, q) for q in range(21, 61) for n in range(max(40, 10 * q), 5001, 7)
                   if not holds(n, q)]
    ok &= not crude_fail
    print(f"  crude lemma (n >= max(40,10q)) failures: {crude_fail[:5]} (count {len(crude_fail)})")

    # exact thresholds n0(q): smallest n >= 4q+1 from which (*) holds for all
    # larger n up to 5000
    print("  exact thresholds n0(q):")
    for q in (2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 50, 80):
        n0 = None
        for n in range(4 * q + 1, 5001):
            if holds(n, q):
                n0 = n
                break
        lat = [n for n in range(n0, 5001) if not holds(n, q)] if n0 else [None]
        ok &= n0 is not None and not lat
        print(f"    q={q:3d}: n0={n0}" + (f"  later failures: {lat[:3]}" if lat else ""))

    # asymptotic threshold gamma = (9 - sqrt(57))/12: verify (*) for all
    # q <= (gamma - 1/50) n in the scan range (exact integer/fraction check)
    from math import sqrt
    gamma = (9 - sqrt(57)) / 12
    g = Fr(101, 1000)  # rational, g < gamma (0.101 < 0.1208)
    asym_fail = []
    for n in range(40, 3001):
        for q in range(2, int(g * n) + 1):
            if not holds(n, q):
                asym_fail.append((n, q))
    ok &= not asym_fail
    print(f"  gamma = {gamma:.6f}; failures for q <= (gamma-1/50)n, n<=3000: {asym_fail[:5]} (count {len(asym_fail)})")

    # max q/n curve
    res = {}
    for n in (100, 500, 2000):
        best = max(q for q in range(2, n // 2) if holds(n, q))
        res[n] = (best, best / n)
    print(f"  max working q at sample n: {res}")

    print("=" * 64)
    print("positive control: q = 2, parity-aware F reproduces r2_density.py")
    print("=" * 64)

    def holds2(n):
        F = n - 1 if n % 2 == 0 else n - 2
        ell = (Fr(1, 2) - Fr(1, n)) / F - Fr(1, 4 * n)
        lhs = ell * (1 - Fr(4, n))
        rhs = (1 + Fr(2, n)) / ((n - 2) * n) + 2 * (1 + Fr(2, n)) / n ** 2
        return lhs > rhs

    odd = [n for n in range(9, 60) if n % 2 == 1 and holds2(n)]
    even = [n for n in range(9, 60) if n % 2 == 0 and holds2(n)]
    ctrl_ok = odd and odd[0] == 19 and even and even[0] == 20 and \
        odd == list(range(19, 60, 2)) and even == list(range(20, 60, 2))
    ok &= ctrl_ok
    print(f"  odd n working from {odd[0] if odd else None}; even from {even[0] if even else None}")

    print("=" * 64)
    print("(S) deletion values: r=2 closed form and sandwich check (exact)")
    print("=" * 64)
    from ml import ml_exact
    bad_cf, bad_sand, bad_eq = [], [], []
    for n in range(4, 31):
        for r in range(1, n):
            V = [i for i in range(1, n) if i != r]
            val, _ = ml_exact(V)
            if val < Fr(1, n - 1):
                bad_sand.append((n, r))
            if val == Fr(1, n - 1) and r != n - 1:
                bad_eq.append((n, r))
            if r == 2 and n >= 5 and val != Fr(2, 2 * (n // 2) + 3):
                bad_cf.append(n)
    ok &= not (bad_cf or bad_sand or bad_eq)
    print(f"  LR([n-1]\\{{r}}) >= 1/(n-1), all (n,r), n in [4,30]: failures {bad_sand[:5]} (count {len(bad_sand)})")
    print(f"  equality 1/(n-1) attained only at r=n-1: violations {bad_eq[:5]} (count {len(bad_eq)})")
    print(f"  r=2 closed form 2/(2*floor(n/2)+3), n in [5,30]: mismatches {bad_cf}")

    print("=" * 64)
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
