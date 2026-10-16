# -*- coding: utf-8 -*-
"""
Independent cross-verification of the paper's computational claims.

The results in the paper were obtained by two routines:
  (A) `mu_criterion.W_intervals` -- builds U(n,r) by intersecting the "good sets"
      of the individual runners (interval-intersection algorithm);
  (A') `ml.ml_exact`             -- evaluates f_V on the breakpoint set T(V).

This file re-derives the same quantities along two DIFFERENT algorithmic paths,
sharing no code with the above (everything is re-implemented from scratch):

  (B) WALL-SUBDIVISION.  Collect every "wall" -- a t with ||i t|| = 1/n for some
      surviving runner i -- sort them, and test the exact midpoint of each
      elementary cell for membership in U(n,r).  Maximal runs of good cells are
      the components of U(n,r).  No set intersection is performed.

  (C) ADAPTIVE BRANCH-AND-BOUND.  Prove LR(V) <= 1/n rigorously WITHOUT
      enumerating candidate maximisers.  On a cell [a,b] the Lipschitz bound
              f_V(t) <= min_i ( ||v_i m|| + v_i (b-a)/2 ),   m = (a+b)/2,
      holds for every t in [a,b].  Cells whose bound is < 1/n + gap are
      discharged, the rest are bisected.  Here
              gap = 1/(n Q),   Q = 2 max(V),
      and the only external input is that LR(V) is a rational with denominator
      at most Q; hence LR(V) != 1/n would force |LR(V) - 1/n| >= 1/(nQ) = gap,
      so LR(V) < 1/n + gap together with f_V(1/n) = 1/n gives LR(V) = 1/n.

Both (B) and (C) use exact rational arithmetic.  A sample of values is
additionally recomputed with sympy.Rational, i.e. a second, independent
rational-number implementation.
"""
from fractions import Fraction as F
from math import gcd
import os
import sys


# ----------------------------------------------------------------- basics ----
def frac_norm(x: F) -> F:
    """||x||, the distance from x to the nearest integer.  Re-implemented."""
    r = x - (x.numerator // x.denominator)      # fractional part in [0,1)
    return r if r + r <= 1 else 1 - r


def f_min(V, t: F) -> F:
    """f_V(t) = min_i ||v_i t||."""
    best = None
    for v in V:
        d = frac_norm(F(v) * t)
        if best is None or d < best:
            best = d
    return best


# ------------------------------------------- (B) wall subdivision of U ------
def walls(n, r):
    """All t in (0,1) with ||i t|| = 1/n for some i in [n-1] \\ {r}."""
    out = set()
    for i in range(1, n):
        if i == r:
            continue
        # i t = b +- 1/n   =>   t = (b n +- 1)/(i n)
        for b in range(0, i + 1):
            for sgn in (1, -1):
                num = b * n + sgn
                if 0 < num < i * n:
                    out.add(F(num, i * n))
    return sorted(out)


def U_components_B(n, r):
    """Components of U(n,r) via wall subdivision + exact midpoint tests."""
    ws = walls(n, r)
    pts = [F(0)] + ws + [F(1)]
    thr = F(1, n)
    good = []
    for j in range(len(pts) - 1):
        a, b = pts[j], pts[j + 1]
        if a == b:
            continue
        mid = (a + b) / 2
        ok = all(frac_norm(F(i) * mid) > thr
                 for i in range(1, n) if i != r)
        good.append((a, b, ok))
    comps, cur = [], None
    for a, b, ok in good:
        if ok:
            cur = (cur[0], b) if cur else (a, b)
        else:
            if cur:
                comps.append(cur)
            cur = None
    if cur:
        comps.append(cur)
    return comps


def F_rep(n, r, x):
    """Largest legal runner  i <= n-1, i != r, i = x (mod r)."""
    reps = [i for i in range(1, n) if i != r and (i - x) % r == 0]
    return max(reps) if reps else None


def I_paper(n, r):
    """I(n,r) as defined in the paper (eq. 1.2), re-implemented."""
    best = None
    for u in range(1, r):
        if gcd(u, r) != 1:
            continue
        for x in (u, r - u):
            v = F_rep(n, r, x)
            if v is None:
                continue
            if best is None or v < best:
                best = v
    return best


def ell_formula(n, r):
    """Longest component length of U(n,r).

    Two regimes.  For 2r <= n-1 the component around p/r is punctured by the
    runner 2r, giving length (2s - I)/(2 r n I).  For 2r > n-1 the component is
    a full interval whose two radii are governed by the classes of u and r-u
    SEPARATELY, so the length is

        (s/(r n)) * max_u ( 1/F(u) + 1/F(r-u) ),

    which is in general NOT 2s/(r n I): the minimum defining I may be attained
    on one side only.
    """
    s, I = n - r, I_paper(n, r)
    if I is None:
        return None
    if 2 * r > n - 1:
        best = None
        for u in range(1, r):
            if gcd(u, r) != 1:
                continue
            a, b = F_rep(n, r, u), F_rep(n, r, r - u)
            if a is None or b is None:
                continue
            cand = F(s, r * n) * (F(1, a) + F(1, b))
            if best is None or cand > best:
                best = cand
        return best
    if 2 * s - I <= 0:
        return F(0)
    return F(2 * s - I, 2 * r * n * I)


def sup_norm_rt(n, r, comps):
    """sup over U(n,r) of ||r t||, computed from the components found by (B)."""
    best = F(0)
    for a, b in comps:
        for t in (a, b):
            d = frac_norm(F(r) * t)
            if d > best:
                best = d
    return best


# --------------------------- (C) branch-and-bound upper bound on LR(V) ------
def prove_LR_le(V, n, max_cells=40_000_000):
    """
    Rigorously certify LR(V) <= 1/n (given f_V(1/n) = 1/n this gives equality),
    without enumerating candidate maximisers.  Returns (True, #cells) or
    (False, offending cell).
    """
    V = sorted(V)
    M = V[-1]
    Q = 2 * M
    gap = F(1, n * Q)
    limit = F(1, n) + gap
    stack = [(F(0), F(1, 2))]          # f_V is symmetric about 1/2
    cells = 0
    while stack:
        a, b = stack.pop()
        cells += 1
        if cells > max_cells:
            return False, ("cell budget exceeded", a, b)
        mid = (a + b) / 2
        half = (b - a) / 2
        # Lipschitz upper bound: every runner gives one, take the smallest
        ub = None
        for v in V:
            cand = frac_norm(F(v) * mid) + F(v) * half
            if ub is None or cand < ub:
                ub = cand
        if ub < limit:
            continue                    # cell discharged
        if half < F(1, 1) / (Q * n * M * 8):
            return False, ("no convergence", a, b, ub)
        stack.append((a, mid))
        stack.append((mid, b))
    return True, cells


# ------------------------------------------------------------- sympy check --
def sympy_spotcheck(V, n, ts):
    """Recompute f_V at given times with sympy.Rational (independent library)."""
    try:
        from sympy import Rational, floor
    except ImportError:
        return None
    out = []
    for t in ts:
        tt = Rational(t.numerator, t.denominator)
        vals = []
        for v in V:
            x = v * tt
            fr = x - floor(x)
            vals.append(min(fr, 1 - fr))
        out.append(min(vals))
    return out


# ------------------------------------------------------------------- main ---
TIGHT_SINGLE = [  # (n, r, w) -- Table 1 of the paper
    (5, 2, 7), (6, 2, 9), (8, 6, 12), (14, 12, 24), (20, 18, 36),
    (26, 24, 48), (32, 30, 60), (32, 30, 90), (33, 30, 60),
    (38, 36, 72), (44, 42, 84),
]
OTHER_TIGHT = [  # further known tight sets
    ((1, 2, 3, 4, 5, 6, 7), 8),
    ((1, 4, 5, 6, 7, 11, 13), 8),
    (tuple(list(range(1, 12)) + [13, 24]), 14),
]


def main():
    print("=" * 74)
    print("(B) wall-subdivision vs (A) interval-intersection vs closed formula")
    print("=" * 74)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mu_criterion import W_intervals          # method (A)

    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    bad_set = bad_len = bad_sup = tested = 0
    for n in range(5, nmax + 1):
        for r in range(2, n):
            compsB = U_components_B(n, r)
            compsA = W_intervals(n, r)
            tested += 1
            if [(a, b) for a, b, in [(x[0], x[1]) for x in compsB]] != \
               [(a, b) for a, b in compsA]:
                bad_set += 1
                if bad_set <= 5:
                    print(f"  COMPONENT MISMATCH n={n} r={r}")
                    print(f"    (A) {compsA}")
                    print(f"    (B) {compsB}")
            lmaxB = max((b - a for a, b in compsB), default=F(0))
            pred = ell_formula(n, r)
            if pred is not None and lmaxB != pred:
                bad_len += 1
                if bad_len <= 8:
                    print(f"  LENGTH MISMATCH n={n} r={r}: (B)={lmaxB} "
                          f"formula={pred}")
            # sup_U ||r t|| = s/(n I)   (used for the Goddyn-Wong margin)
            if 2 * r > n - 1 and compsB:
                I = I_paper(n, r)
                if I is not None:
                    if sup_norm_rt(n, r, compsB) != F(n - r, n * I):
                        bad_sup += 1
                        if bad_sup <= 8:
                            print(f"  SUP MISMATCH n={n} r={r}: "
                                  f"(B)={sup_norm_rt(n, r, compsB)} "
                                  f"formula={F(n - r, n * I)}")
    print(f"  pairs tested: {tested}   component-set mismatches: {bad_set}   "
          f"length mismatches: {bad_len}   sup mismatches: {bad_sup}")

    print()
    print("=" * 74)
    print("(C) branch-and-bound certification of LR(V) = 1/n (no breakpoints)")
    print("=" * 74)
    allV = [(tuple(sorted([x for x in range(1, n) if x != r] + [w])), n)
            for (n, r, w) in TIGHT_SINGLE] + OTHER_TIGHT
    for V, n in allV:
        at = f_min(V, F(1, n))
        ok, info = prove_LR_le(V, n)
        tag = "CERTIFIED LR = 1/n" if (ok and at == F(1, n)) else "FAILED"
        print(f"  n={n:3d} max={max(V):3d}  f(1/n)={at}  cells={info}  -> {tag}")

    print()
    print("=" * 74)
    print("(D) sympy.Rational spot-check (independent rational implementation)")
    print("=" * 74)
    for V, n in allV[:6]:
        ts = [F(1, n), F(2, n), F(1, 2), F(1, max(V))]
        mine = [f_min(V, t) for t in ts]
        theirs = sympy_spotcheck(V, n, ts)
        if theirs is None:
            print("  sympy unavailable")
            break
        agree = all(F(str(a)) == F(int(b.p), int(b.q))
                    for a, b in zip(mine, theirs))
        print(f"  n={n:3d}  values {[str(x) for x in mine]}  agree={agree}")


if __name__ == "__main__":
    main()
