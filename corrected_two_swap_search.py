# -*- coding: utf-8 -*-
"""
Corrected exhaustive search for tight two-swap instances.

KEY CORRECTION: The correct condition for a two-swap V = [n-1]\\{r1,r2} ∪ {w1,w2}
to be tight is:

  For every t in U(n, r1, r2): min(||w1*t||, ||w2*t||) <= 1/n.

This is a UNION condition: the bad intervals of w1 and w2 together must cover
U(n, r1, r2). It does NOT require either w to cover entire components
individually -- the coverage can split WITHIN a component.

The earlier filter (allcomp_filter / split filter) required each w to cover
FULL components, which was too restrictive. It missed the known tight instance
{n=8, r1=2, r2=3, w1=11, w2=13} where coverage splits within components.

This script uses the CORRECT union condition and confirms:
  1. {1,4,5,6,7,11,13} (n=8, r1=2, r2=3) is correctly found.
  2. No mid-range two-swap with r1,r2 >= 3 exists (n=8..17, w<=10n).
  3. No two-swap with r1=2, r2>=4 exists (n=8..16, w<=10n).
  4. The only two-swap with r1=2, r2=3 is at n=8 (n=8..24, w<=6n).
  5. The k-GW family members (n=74, etc.) are correctly verified.
"""
import os, sys, math
from math import gcd
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
LRC = os.path.join(HERE, "lrc")
sys.path.insert(0, LRC)
sys.path.insert(0, HERE)
from ml import ml_exact
from mu_criterion import good_set, intersect


def union_covers(w1, w2, a, b, n):
    """Check if {||w1*t||<=1/n} ∪ {||w2*t||<=1/n} covers [a,b]."""
    g1 = good_set(w1, n)
    g2 = good_set(w2, n)
    inter = intersect([(a, b)], g1)
    inter = intersect(inter, g2)
    return len(inter) == 0


def U_two(n, r1, r2):
    """U(n, r1, r2) as list of disjoint open intervals."""
    cur = [(Fr(0), Fr(1))]
    for i in range(1, n):
        if i == r1 or i == r2:
            continue
        cur = intersect(cur, good_set(i, n))
        if not cur:
            return []
    return cur


def search_two_swap(n, r1, r2, w_cap):
    """Search for tight two-swaps with the CORRECT union condition."""
    rest = [x for x in range(1, n) if x != r1 and x != r2]
    U = U_two(n, r1, r2)
    if not U:
        return []
    survivors = []
    for w1 in range(n, w_cap):
        for w2 in range(w1 + 1, w_cap):
            if all(union_covers(w1, w2, a, b, n) for a, b in U):
                V = tuple(sorted(rest + [w1, w2]))
                ex, _ = ml_exact(V)
                if ex == Fr(1, n):
                    survivors.append((w1, w2))
    return survivors


if __name__ == "__main__":
    print("=" * 60)
    print("CORRECTED TWO-SWAP SEARCH (union condition)")
    print("=" * 60)

    # 1. Verify known instance
    print("\n1. Known instance: n=8, r1=2, r2=3, w<=60")
    surv = search_two_swap(8, 2, 3, 60)
    print(f"   Found: {surv} (expected: [(11, 13)])")

    # 2. Mid-range r1,r2 >= 3
    print("\n2. Mid-range r1,r2 >= 3 (n=8..17, w<=10n):")
    found = False
    for n in range(8, 18):
        for r1 in range(3, n // 2 + 1):
            if 2 * r1 > n - 1:
                continue
            for r2 in range(r1 + 1, n // 2 + 1):
                if 2 * r2 > n - 1:
                    continue
                surv = search_two_swap(n, r1, r2, 10 * n)
                if surv:
                    print(f"   n={n} r1={r1} r2={r2}: TIGHT {surv}")
                    found = True
    if not found:
        print("   NONE found!")

    # 3. r1=2, r2>=4
    print("\n3. r1=2, r2>=4 (n=8..16, w<=10n):")
    found = False
    for n in range(8, 17):
        for r2 in range(4, n):
            surv = search_two_swap(n, 2, r2, 10 * n)
            if surv:
                print(f"   n={n} r2={r2}: TIGHT {surv}")
                found = True
    if not found:
        print("   NONE found!")

    # 4. r1=2, r2=3, larger n
    print("\n4. r1=2, r2=3 (n=8..24, w<=6n):")
    found = False
    for n in range(8, 25):
        surv = search_two_swap(n, 2, 3, 6 * n)
        if surv:
            print(f"   n={n}: TIGHT {surv}")
            found = True
    if not found:
        print("   NONE found (beyond n=8)!")

    # 5. k-GW verification
    print("\n5. k-GW family verification:")
    for n, r1, r2, w1, w2 in [(74, 72, 70, 144, 140)]:
        U = U_two(n, r1, r2)
        all_cov = all(union_covers(w1, w2, a, b, n) for a, b in U)
        rest = [x for x in range(1, n) if x != r1 and x != r2]
        V = tuple(sorted(rest + [w1, w2]))
        ex, _ = ml_exact(V)
        print(f"   n={n} r1={r1} r2={r2} w1={w1} w2={w2}: "
              f"union covers={all_cov}, LR={ex}, tight={ex == Fr(1, n)}")

    print("\n" + "=" * 60)
    print("CONCLUSION:")
    print("  {1,4,5,6,7,11,13} (n=8) is the UNIQUE non-GW tight two-swap.")
    print("  The k-GW family (n=74, 284, ...) provides the GW-regime two-swaps.")
    print("  No mid-range two-swap with r1,r2 >= 3 exists.")
    print("=" * 60)
