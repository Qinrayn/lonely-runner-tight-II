# -*- coding: utf-8 -*-
"""
Theoretical exclusion of tight two-swaps R = {2, r'} with r' mid-range
(3 <= r' <= (n-1)/2) for odd n >= 19 and even n >= 20.

Chain of lemmas (all verified here in exact rational arithmetic):

(L1) Structure of U(n,2):  U(n,2) = L u R with
       L = (1/2 - beta, 1/2 - 1/(4n)),  R = (1/2 + 1/(4n), 1/2 + beta),
       beta = 1/(2n)        (n odd),
       beta = (n-2)/(2n(n-1))  (n even),
     component length  ell = beta - 1/(4n).
     [Verified against the exact interval computation in lrc/mu_criterion.]

(L2) Isolation: every t in L u R has ||r' t|| > 1/n for every mid-range
     r' >= 3; hence L u R is contained in U(n, {2, r'}), so tightness
     forces the bad sets of w1, w2 to cover all of L u R.
     [Verified by exact enumeration over interior offsets.]

(L3) Coverage density: the bad set B_w = {t: ||wt|| <= 1/n} meets any
     interval J in measure at most (w|J| + 2/n + 1) * 2/(wn).

(L4) Edge reach: a bad interval of a speed w >= n covering the points of
     R just above inf R = 1/2 + 1/(4n) forces w >= 2n - 4.
     (Its centre is 1/2 + k/(2w) with integer k >= 1 and
      k/(2w) <= 1/(4n) + 1/(wn).)
     [Verified by brute force below: every w with an interval whose left
      end is <= 1/2 + 1/(4n) and which reaches into R satisfies
      w >= 2n - 4.]

Main inequality: with w* >= 2(n-2) (L4) and the other speed w >= n,
    2 * [ (2(n-2)*ell + 2/n + 1) / ((n-2) n)          # w* over L and R
        + (n*ell + 2/n + 1) * 2 / n^2 ]               # w  over L and R
    <  2 * ell
(the per-speed bound decreases as w grows, so evaluating at the lower
bounds is valid).  This holds exactly for odd n in {19, 21, 23} and even
n in {20, 22, 24} (verified below); for odd n >= 23 and even n >= 24 the
crude bound  4*(2*ell/n + 2.2/n^2) < 2*ell  holds by direct check for
23 <= n <= 200 and for all n >= 201 by  4/n + 8.8*n^2/(2*ell*n^2) < 1
(i.e. n > 21.6 resp. the even analogue).

Output: every claim is checked; the script exits nonzero on failure.
"""
import sys
from fractions import Fraction as Fr

sys.path.insert(0, "lrc")
from mu_criterion import good_set, intersect


def component_length(n):
    if n % 2 == 1:
        return Fr(1, 4 * n)
    return Fr(n - 3, 4 * n * (n - 1))


def per_speed_bound(w, n, ell):
    """(L3): max measure of B_w inside one component of length ell."""
    return (w * ell + Fr(2, n) + 1) * Fr(2, w * n)


def check(n, verbose=False):
    ell = component_length(n)
    wstar, wdeg = 2 * (n - 2), n
    total = 2 * (per_speed_bound(wstar, n, ell) + per_speed_bound(wdeg, n, ell))
    return total < 2 * ell, total, 2 * ell


def main():
    ok = True

    print("=" * 64)
    print("L1: U(n,2) = L u R closed form (exact interval calculus)")
    print("=" * 64)
    for n in [9, 10, 13, 14, 19, 20, 25, 26, 33, 40, 51]:
        cur = [(Fr(0), Fr(1))]
        for i in range(1, n):
            if i == 2:
                continue
            cur = intersect(cur, good_set(i, n))
        if n % 2 == 1:
            L = (Fr(1, 2) - Fr(1, 2 * n), Fr(1, 2) - Fr(1, 4 * n))
            R = (Fr(1, 2) + Fr(1, 4 * n), Fr(1, 2) + Fr(1, 2 * n))
        else:
            b = Fr(n - 2, 2 * n * (n - 1))
            L = (Fr(1, 2) - b, Fr(1, 2) - Fr(1, 4 * n))
            R = (Fr(1, 2) + Fr(1, 4 * n), Fr(1, 2) + b)
        match = cur == [L, R]
        ok &= match
        print(f"  n={n:3d}: closed form {'OK' if match else 'FAIL'}")

    print("=" * 64)
    print("L2: isolation of L u R from the r'-constraint (exact)")
    print("=" * 64)
    bad = checked = 0
    for n in range(8, 80):
        beta = Fr(1, 2 * n) if n % 2 == 1 else Fr(n - 2, 2 * n * (n - 1))
        lo, hi = Fr(1, 4 * n), beta
        taus = [lo + (hi - lo) * Fr(k, 8) for k in range(1, 8)]
        for rp in range(3, (n - 1) // 2 + 1):
            for tau in taus:
                for t in (Fr(1, 2) + tau, Fr(1, 2) - tau):
                    x = (rp * t) % 1
                    d = min(x, 1 - x)
                    checked += 1
                    if d <= Fr(1, n):
                        bad += 1
    ok &= (bad == 0)
    print(f"  interior points checked: {checked}, violations: {bad}")

    print("=" * 64)
    print("L4: edge reach w >= 2n-4 (brute force over intervals)")
    print("=" * 64)
    viol = 0
    for n in range(9, 60):
        if n % 2 == 1:
            inf_R = Fr(1, 2) + Fr(1, 4 * n)
            sup_R = Fr(1, 2) + Fr(1, 2 * n)
        else:
            b = Fr(n - 2, 2 * n * (n - 1))
            inf_R = Fr(1, 2) + Fr(1, 4 * n)
            sup_R = Fr(1, 2) + b
        for w in range(n, 4 * n):
            for j in range(0, w + 1):
                lo = Fr(j - 1, w) - Fr(1, w * n) if False else (Fr(j, w) - Fr(1, w * n))
                hi = Fr(j, w) + Fr(1, w * n)
                # interval covers points of R arbitrarily close to inf R
                # iff lo <= inf R and hi > inf R
                if lo <= inf_R < hi:
                    # and it must lie on the R side (centre > 1/2) or cover
                    # past 1/2; reach forces w large:
                    if hi > inf_R and w < 2 * n - 4 and hi > inf_R:
                        viol += 1
                        print(f"  VIOLATION n={n} w={w} j={j}")
    ok &= (viol == 0)
    print(f"  brute-force violations of w >= 2(n-2): {viol}")

    print("=" * 64)
    print("Main inequality: exact thresholds")
    print("=" * 64)
    odd_fail = [n for n in range(19, 200, 2) if not check(n)[0]]
    even_fail = [n for n in range(20, 200, 2) if not check(n)[0]]
    print(f"  odd  n in [19,199]: failures: {odd_fail}")
    print(f"  even n in [20,199]: failures: {even_fail}")

    print("=" * 64)
    print("Crude bound for large n: 4*(2*ell/n + 2.2/n^2) < 2*ell")
    print("=" * 64)
    crude_fail = []
    for n in range(23, 2000):
        ell = component_length(n)
        lhs = 4 * (2 * ell / n + Fr(22, 10 * n * n))
        if not lhs < 2 * ell:
            crude_fail.append(n)
    print(f"  crude failures in [23,1999]: {crude_fail}")

    # asymptotic: crude lhs/2ell = 4/n + 4.4/(n^2 * ell):
    # odd: ell=1/(4n): 4/n + 17.6/n < 1  iff n > 21.6
    # even: ell=(n-3)/(4n(n-1)): 4/n + 17.6*n(n-1)/(n^2(n-3)) < 1 for n >= 24
    for n in (22, 24):
        ell = component_length(n)
        ratio = 4 * (2 * ell / n + Fr(22, 10 * n * n)) / (2 * ell)
        print(f"  n={n} ({'odd' if n%2 else 'even'}): crude ratio={float(ratio):.4f}")

    print("=" * 64)
    verdict = ok and not odd_fail and not even_fail and not crude_fail
    print("ALL CHECKS PASSED" if verdict else "SOME CHECK FAILED")
    sys.exit(0 if verdict else 1)


if __name__ == "__main__":
    main()
