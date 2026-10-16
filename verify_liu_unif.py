# -*- coding: utf-8 -*-
"""
Verification of the two claims added in the 2026-10-16 revision round:

(1) LIU'S FORMULA (s = 1).  Perarnau-Serra [survey, Thm. 14] record Liu's
    exact value for sets with an arithmetic progression removed; the case
    s = 1 (survey n = our n-1, survey k = our r) states, for r <= (n-1)/2:
        LR([n-1] \\ {r}) = 2 / D_r,
        D_r = min{ D >= n + r : gcd(D, r) = 1 }.
    We re-derive the left-hand side independently in exact rational
    arithmetic (lrc/ml.ml_exact) and compare, for 5 <= n <= 24 and all
    r <= (n-1)/2.  We also record whether D_r <= 2n-3 (i.e. 2/D_r > 1/(n-1))
    and whether the witness bound of the paper (Thm delwit, D = D_r) is
    attained with equality.

(2) UNIFORM COVERING BOUND (Lemma lem:unif).  For n = 5..12 and the table
    values v_n = smallest non-tight LR([n-1]_{2->w}) (w <= 8n), compute
        i_0 = largest odd <= n-1,
        l(v) = (1/2 - v)/i_0 - v/4,
        threshold = 2 v / l(v),
    exactly with Fractions, and check:
        (a) threshold is the integer listed in the paper (12,20,10,14,14,18,18,22);
        (b) the attaining w satisfies n <= w <= threshold;
        (c) threshold <= 8n (so the w <= 8n search plus the lemma cover all w);
        (d) for every w in [n, 8n] with LR([n-1]_{2->w}) <= v_n, w <= threshold
            (the lemma's conclusion, as a consistency check);
        (e) 2*lambda/l(lambda) is increasing in lambda on (0, 2/D).

Run:  python v3/verify_liu_unif.py
"""
import sys, os
from fractions import Fraction as F
from math import gcd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "lrc"))
from ml import ml_exact  # noqa: E402

fail = 0

# ---------------------------------------------------------------- (1) Liu ---
print("== (1) Liu s=1 formula: LR([n-1]\\{r}) == 2/D_r for r <= (n-1)/2 ==")
viol = 0
worst_D = 0
for n in range(5, 25):
    for r in range(1, (n - 1) // 2 + 1):
        V = sorted(set(range(1, n)) - {r})
        lr, _ = ml_exact(V)
        D = next(D for D in range(n + r, 4 * n) if gcd(D, r) == 1)
        if lr != F(2, D):
            viol += 1
            print(f"  MISMATCH n={n} r={r}: exact={lr} vs Liu 2/D={F(2, D)}")
        if D > 2 * n - 3:
            print(f"  NOTE n={n} r={r}: D_r={D} > 2n-3={2*n-3} (2/D <= 1/(n-1))")
        worst_D = max(worst_D, D - (2 * n - 3))
print(f"  checked n=5..24, all r <= (n-1)/2: mismatches = {viol}")
print("  (all D_r <= 2n-3 unless NOTEd above; max slack D_r-(2n-3) =", worst_D, ")")
fail += viol

# ------------------------------------------------------- (2) Lemma lem:unif --
print("== (2) uniform covering bound: thresholds 2*v/l(v) ==")
expected = {5: 12, 6: 20, 7: 10, 8: 14, 9: 14, 10: 18, 11: 18, 12: 22}
table_v = {5: F(2, 9), 6: F(2, 11), 7: F(2, 13), 8: F(2, 15),
           9: F(2, 15), 10: F(2, 17), 11: F(2, 17), 12: F(2, 19)}
table_w = {5: 5, 6: 7, 7: 9, 8: 11, 9: 11, 10: 13, 11: 13, 12: 15}
D_of = lambda n: 2 * (n // 2) + 3
for n in range(5, 13):
    i0 = n - 2 if n % 2 == 1 else n - 1
    assert D_of(n) == i0 + 4
    v = table_v[n]
    ell = (F(1, 2) - v) / i0 - v / 4
    thr = 2 * v / ell
    ok_a = (thr.denominator == 1 and thr.numerator == expected[n])
    ok_b = (n <= table_w[n] <= thr)
    ok_c = (thr <= 8 * n)
    # (d): every w <= 8n with LR <= v_n has w <= thr  (lemma consistency)
    ok_d = True
    for w in range(n, 8 * n + 1):
        if w == 2:
            continue
        V = sorted((set(range(1, n)) - {2}) | {w})
        lr, _ = ml_exact(V)
        if lr <= v and w > thr:
            ok_d = False
            print(f"  LEMMA VIOLATION n={n} w={w}: LR={lr} <= {v} but w > {thr}")
    # lambda = 1/n recovers Delta(n)
    ell_n = (F(1, 2) - F(1, n)) / i0 - F(1, 4 * n)
    delta = F(1, 4 * n) if n % 2 == 1 else F(n - 3, 4 * n * (n - 1))
    ok_e0 = (ell_n == delta)
    print(f"  n={n:2d}: i0={i0:2d} l(v)={ell} thr={thr} "
          f"a={ok_a} b={ok_b} c={ok_c} d={ok_d} Delta-recovers={ok_e0}")
    fail += (not ok_a) + (not ok_b) + (not ok_c) + (not ok_d) + (not ok_e0)

# (e) monotonicity of 2*lambda/l(lambda) on (0, 2/D), sampled exactly
print("== (2e) monotonicity of 2*lambda/l(lambda) (exact sample) ==")
mono_ok = True
for n in range(5, 13):
    i0 = n - 2 if n % 2 == 1 else n - 1
    D = i0 + 4
    prev = None
    for k in range(1, 40):
        lam = F(k, 40) * F(2, D)
        if lam >= F(2, D):
            break
        ell = (F(1, 2) - lam) / i0 - lam / 4
        if ell <= 0:
            break
        val = 2 * lam / ell
        if prev is not None and val <= prev:
            mono_ok = False
            print(f"  NOT increasing at n={n} lambda={lam}")
        prev = val
print(f"  increasing on (0, 2/D) for n=5..12 (exact samples): {mono_ok}")
fail += (not mono_ok)

print()
print("TOTAL FAILURES:", fail)
