# -*- coding: utf-8 -*-
"""
The mid-range classification as a finite check.

By Theorem 1.6 of paper 1 (effective bound), a tight [n-1]_{r->w} with
2r <= n-1 forces  n <= 6r  and  w <= 4 r I / (2s - I),  s = n-r,
I = min over units u mod r of min(F(u), F(r-u)).

Hence for each fixed r the whole classification is the finite scan
    n in [2r+1, 6r],   w in [n, floor(4rI/(2s-I))],
and the conjecture (Theorem A) is that for r >= 3 it finds nothing,
while r = 2 yields exactly (5,7) and (6,9).

Exactness: float prescreen on the breakpoint set (provably safe, Lemma 5.1 of
paper 1), exact rational confirmation for survivors.
"""
import os, sys, time, json
from math import gcd
from fractions import Fraction as F
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
LRC = os.path.join(HERE, "lrc")
sys.path.insert(0, LRC)
from ml import ml_exact                      # exact LR via breakpoints
from twoswap_hunt import certifies_not_tight  # numba float prescreen


def Frep(n, r, x):
    reps = [i for i in range(1, n) if i != r and (i - x) % r == 0]
    return max(reps) if reps else None


def Ival(n, r):
    best = None
    for u in range(1, r):
        if gcd(u, r) != 1:
            continue
        for x in (u, r - u):
            v = Frep(n, r, x)
            if v is not None and (best is None or v < best):
                best = v
    return best


def classify_r(r, log=print):
    """Complete finite classification of tight [n-1]_{r->w} with 2r <= n-1."""
    hits, checked = [], 0
    t0 = time.time()
    for n in range(2 * r + 1, 6 * r + 1):
        s = n - r
        I = Ival(n, r)
        if I is None or 2 * s - I <= 0:
            continue                      # U empty is impossible (GW Lem 2.2) but guard
        wmax = (4 * r * I) // (2 * s - I)
        if wmax < n:
            continue
        rest = [x for x in range(1, n) if x != r]
        thr = 1.0 / n + 1e-9
        for w in range(n, wmax + 1):
            V = sorted(rest + [w])
            checked += 1
            Va = np.array(V, dtype=np.int64)
            if certifies_not_tight(Va, n - 1, thr):
                continue
            ex, _ = ml_exact(tuple(V))
            if ex == F(1, n):
                hits.append((n, w))
                log(f"    TIGHT: r={r} n={n} w={w}")
            elif ex < F(1, n):
                log(f"    !!! LRC COUNTEREXAMPLE r={r} n={n} w={w} LR={ex}")
    dt = time.time() - t0
    log(f"  r={r:3d}: n in [{2*r+1},{6*r}], {checked} candidates, "
        f"{len(hits)} tight, {dt:.0f}s")
    return {"r": r, "checked": checked, "tight": hits, "seconds": round(dt, 1)}


if __name__ == "__main__":
    rs = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 \
        else list(range(2, 13))
    out = []
    for r in rs:
        out.append(classify_r(r))
    path = os.path.join(HERE, f"midrange_r{rs[0]}_{rs[-1]}.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=1)
    print("\nSUMMARY")
    for rec in out:
        print(f"  r={rec['r']:3d}: tight = {rec['tight'] if rec['tight'] else 'NONE'}"
              f"  ({rec['checked']} candidates)")
