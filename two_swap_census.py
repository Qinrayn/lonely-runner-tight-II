# -*- coding: utf-8 -*-
"""
Exhaustive census of tight two-swap modifications
    V = [n-1] \\ {r1, r2}  U  {w1, w2},   n <= w1 < w2 <= 8n,
over ALL removal pairs 1 <= r1 < r2 <= n-1 (any regime), n = 8..N_MAX.

Method (provably-safe screening + exact certification):

  * Components of U(n,r) are computed in floats with the bad intervals
    EXPANDED by EPS_BAD, so the float component list under-approximates
    the true uncovered region: every float component is a genuine
    component-piece of the exact U, hence every float-gap argument is
    sound (points of the float region must be covered in any tight set).
  * A speed w can cover the points of a component C=(a,b) arbitrarily
    close to inf C only if some integer j lies in [w*a - 1/n, w*a + 1/n]
    (left-edge set S_L); similarly the right edge gives S_R.  Both are
    necessary conditions, computed with generous slack (superset).
  * Candidate pairs surviving the edge-set masks get a float coverage
    test; a pair is DISCARDED there only if some component shows a gap
    > GAP_TOL = 1e-9.  True gaps inside a component, when they exist,
    are at least 1/(n * w1 * w2) >= 1/(64 n^3) ~ 5.6e-9 (n=140), far
    above both the float noise (~1e-14) and GAP_TOL, so no genuine
    covering pair is ever discarded.
  * Every survivor is re-verified in EXACT rational arithmetic:
    exact U(n,{r1,r2}) from lrc.mu_criterion, exact interval coverage,
    and LR(V) = 1/n from lrc.ml.ml_exact (breakpoint certification).

Positive controls built in: n=8, R={2,3} must rediscover {11,13};
n=74, R={70,72} must rediscover the k-GW pair {140,144}.
"""
import os
import sys
import time
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "lrc"))
sys.path.insert(0, HERE)
from ml import ml_exact                      # noqa: E402
from mu_criterion import good_set, intersect  # noqa: E402

EPS_BAD = 1e-12   # expansion of bad intervals (float region subset of true)
EPS_SET = 1e-7    # slack in edge-set windows (superset of exact sets)
GAP_TOL = 1e-9    # discard only on gaps strictly larger than this


# ---------------------------------------------------------------- floats --

def uncovered_float(n, r):
    """Components of U(n,r) in floats (under-approximation of the truth)."""
    los, his = [], []
    for i in range(1, n):
        if i == r:
            continue
        j = np.arange(i + 1)
        los.append((j - 1.0 / n) / i - EPS_BAD)
        his.append((j + 1.0 / n) / i + EPS_BAD)
    lo = np.concatenate(los)
    hi = np.concatenate(his)
    order = np.argsort(lo, kind="stable")
    lo, hi = lo[order], hi[order]
    # merge: group boundaries where lo exceeds every earlier hi
    run_hi = np.maximum.accumulate(hi)
    new_bad = np.empty(len(lo), dtype=bool)
    new_bad[0] = True
    new_bad[1:] = lo[1:] > run_hi[:-1]
    starts = np.nonzero(new_bad)[0]
    ends = np.r_[starts[1:] - 1, len(lo) - 1]
    mlo, mhi = lo[starts], run_hi[ends]
    # complement in [0,1]
    comps = []
    cur = 0.0
    for a, b in zip(mlo, mhi):
        a = max(a, 0.0)
        if a > cur + 1e-15:
            comps.append((cur, min(a, 1.0)))
        cur = max(cur, min(b, 1.0))
    if cur < 1.0 - 1e-15:
        comps.append((cur, 1.0))
    return [(a, b) for a, b in comps if b - a > 1e-13]


def union_floats(u1, u2):
    """U(n,{r1,r2}) contains U(n,r1) U U(n,r2) always; when r1+r2 < n the
    bad sets of r1 and r2 are disjoint (|j/r1 - j'/r2| >= 1/(r1 r2) >
    (r1+r2)/(r1 r2 n)) and equality holds.  The union is used as a sound
    under-approximation: pairs failing to cover it cannot cover the true
    region, and every survivor is re-checked against the exact U(n,{r1,r2})."""
    raw = sorted(u1 + u2)
    out = [list(raw[0])]
    for a, b in raw[1:]:
        if a <= out[-1][1] + 1e-13:
            out[-1][1] = max(out[-1][1], b)
        else:
            out.append([a, b])
    return [(a, b) for a, b in out]


# ------------------------------------------------------------- screening --

class TripleData:
    """Edge-reach candidate sets for one (n, r1, r2)."""

    def __init__(self, n, comps, w_lo, w_hi):
        self.n = n
        self.comps = comps
        self.W = np.arange(w_lo, w_hi + 1)
        self.nw = len(self.W)
        inv = 1.0 / n
        self.SL = []   # per component: boolean array over W (left edge)
        self.SR = []
        self.SLs = []  # same as frozenset of indices (fast membership)
        self.SRs = []
        for a, b in comps:
            x = self.W * a
            j = np.ceil(x - inv - EPS_SET)
            sl = (j <= x + inv + EPS_SET) & (j >= 0)
            y = self.W * b
            j2 = np.ceil(y - inv - EPS_SET)
            sr = (j2 <= y + inv + EPS_SET) & (j2 >= 0)
            self.SL.append(sl); self.SR.append(sr)
            self.SLs.append(frozenset(np.nonzero(sl)[0].tolist()))
            self.SRs.append(frozenset(np.nonzero(sr)[0].tolist()))


def float_covers(comp_list, w1, w2, n):
    """True if the float intervals of w1, w2 appear to cover every
    component (gaps <= GAP_TOL tolerated); False if some gap > GAP_TOL."""
    inv = 1.0 / n
    for a, b in comp_list:
        segs = []
        for w in (w1, w2):
            jlo = int(np.ceil(w * a - inv - EPS_SET))
            jhi = int(np.floor(w * b + inv + EPS_SET))
            for jj in range(max(jlo, 0), min(jhi, w) + 1):
                lo = (jj - inv) / w
                hi = (jj + inv) / w
                lo, hi = max(lo, a), min(hi, b)
                if hi > lo:
                    segs.append((lo, hi))
        if not segs:
            return False
        segs.sort()
        cur = a
        for lo, hi in segs:
            if lo > cur + GAP_TOL:
                return False
            cur = max(cur, hi)
        if cur < b - GAP_TOL:
            return False
    return True


# ----------------------------------------------------------- exact stage --

def exact_check(n, r1, r2, w1, w2):
    """Exact verification: does {w1,w2} cover U(n,{r1,r2})?  Returns
    'tight', 'covers' (but not tight), or 'gap'."""
    cur = [(Fr(0), Fr(1))]
    for i in range(1, n):
        if i in (r1, r2):
            continue
        cur = intersect(cur, good_set(i, n))
        if not cur:
            return "empty"
    for a, b in cur:
        segs = []
        for w in (w1, w2):
            for j in range(0, w + 1):
                lo, hi = Fr(j, w) - Fr(1, w * n), Fr(j, w) + Fr(1, w * n)
                lo, hi = max(lo, a), min(hi, b)
                if hi > lo:
                    segs.append((lo, hi))
        segs.sort()
        cov = a
        for lo, hi in segs:
            if lo > cov:
                return "gap"
            cov = max(cov, hi)
        if cov < b:
            return "gap"
    rest = [x for x in range(1, n) if x != r1 and x != r2]
    V = tuple(sorted(rest + [w1, w2]))
    ex, _ = ml_exact(V)
    return "tight" if ex == Fr(1, n) else "covers"


# ------------------------------------------------------------------ main --

def run_census(n_min, n_max, w_mult=8, r_filter=None, deep=False, log=print):
    stats = dict(triples=0, skip_empty=0, skip_edge=0, pairs=0,
                 float_pass=0, exact=0, tight=[])
    for n in range(n_min, n_max + 1):
        t0 = time.time()
        U = {}
        for r in range(1, n):
            U[r] = uncovered_float(n, r)
        W_LO, W_HI = n, w_mult * n
        n_triples = n_pairs = n_exact = 0
        for r1 in range(1, n):
            for r2 in range(r1 + 1, n):
                if r_filter and not r_filter(n, r1, r2):
                    continue
                comps = union_floats(U[r1], U[r2])
                if not comps:
                    stats["skip_empty"] += 1
                    continue
                td = TripleData(n, comps, W_LO, W_HI)
                # any component with an empty edge set kills the triple
                if not all(sl.any() and sr.any() for sl, sr in zip(td.SL, td.SR)):
                    stats["skip_edge"] += 1
                    continue
                stats["triples"] += 1
                n_triples += 1
                # most restrictive components first (fewest edge candidates)
                sizes = [(int(sl.sum()) + int(sr.sum()), c)
                         for c, (sl, sr) in enumerate(zip(td.SL, td.SR))]
                sizes.sort()
                order = [c for _, c in sizes]
                found = census_triple(td, order)
                for w1, w2 in found:
                    n_pairs += 1
                    stats["pairs"] += 1
                    if float_covers(comps, w1, w2, n):
                        stats["float_pass"] += 1
                        n_exact += 1
                        stats["exact"] += 1
                        res = exact_check(n, r1, r2, w1, w2)
                        if res == "tight":
                            stats["tight"].append((n, r1, r2, w1, w2))
                            log(f"  TIGHT: n={n} R={{{r1},{r2}}} "
                                f"W={{{w1},{w2}}}")
                        elif res == "covers":
                            log(f"  note: covers-but-not-tight n={n} "
                                f"R={{{r1},{r2}}} W={{{w1},{w2}}}")
        log(f"n={n}: triples={n_triples} pairs={n_pairs} exact={n_exact} "
            f"({time.time()-t0:.1f}s)")
    return stats


def census_triple(td, order):
    """Candidate pairs passing the edge-set masks on the most restrictive
    components (necessary conditions), lazily verified on all components."""
    out = set()
    W = td.W
    # prefilter component pool
    pool = order[:12]
    c0 = order[0]
    for ia in np.nonzero(td.SL[c0])[0]:
        wa = int(W[ia])
        cand = np.ones(td.nw, dtype=bool)
        cand[ia] = False
        ok = True
        for c in pool:
            sla = bool(td.SL[c][ia])
            sra = bool(td.SR[c][ia])
            cand &= td.SL[c] | sla
            cand &= td.SR[c] | sra
            if not cand.any():
                ok = False
                break
        if not ok:
            continue
        for ib in np.nonzero(cand)[0]:
            wb = int(W[ib])
            w1, w2 = min(wa, wb), max(wa, wb)
            if (w1, w2) in out:
                continue
            if full_mask_ok(td, ia, ib):
                out.add((w1, w2))
    # pairs where the C0-left-edge coverer is w_b: symmetric, covered by
    # iterating wa over SL[c0] and partners over all (done above).
    return out


def full_mask_ok(td, ia, ib):
    for sl, sr in zip(td.SLs, td.SRs):
        if ia not in sl and ib not in sl:
            return False
        if ia not in sr and ib not in sr:
            return False
    return True


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--nmax", type=int, default=140)
    ap.add_argument("--nmin", type=int, default=8)
    ap.add_argument("--wmult", type=int, default=8)
    ap.add_argument("--deep", action="store_true",
                    help="deep window: R={2,r'} for 9<=n<=18, w<=40n")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()

    logfile = open("two_swap_census.log", "a", encoding="utf-8")

    def log(msg):
        print(msg)
        logfile.write(msg + "\n")
        logfile.flush()

    if args.smoke:
        log("=== SMOKE TEST n=8..12 (expect TIGHT n=8 R={2,3} W={11,13}) ===")
        st = run_census(8, 12, log=log)
        log(f"smoke summary: {st}")
        assert any(t[0] == 8 and t[1] == 2 and t[2] == 3 for t in st["tight"]), \
            "positive control n=8 {11,13} NOT found"
        log("smoke OK")
        return

    if args.deep:
        log("=== DEEP WINDOW: R={2,r'}, n=9..18, w<=40n ===")
        st = run_census(9, 18, w_mult=40,
                        r_filter=lambda n, r1, r2:
                        (r1 == 2 or r2 == 2), log=log)
        log(f"deep summary: triples={st['triples']} pairs={st['pairs']} "
            f"exact={st['exact']} tight={st['tight']}")
        return

    log("=== POSITIVE CONTROL: n=74 k-GW (expect TIGHT R={70,72} W={140,144}) ===")
    st = run_census(74, 74, r_filter=lambda n, r1, r2: (r1, r2) == (70, 72),
                    log=log)
    assert any(t[0] == 74 for t in st["tight"]), \
        "positive control n=74 {140,144} NOT found"
    log("control OK")

    log(f"=== FULL CENSUS: all pairs, n={args.nmin}..{args.nmax}, "
        f"w<={args.wmult}n ===")
    t0 = time.time()
    st = run_census(args.nmin, args.nmax, w_mult=args.wmult, log=log)
    log(f"TOTAL ({time.time()-t0:.0f}s): triples={st['triples']} "
        f"skip_empty={st['skip_empty']} skip_edge={st['skip_edge']} "
        f"pairs={st['pairs']} float_pass={st['float_pass']} "
        f"exact={st['exact']}")
    log(f"TIGHT SETS FOUND: {st['tight']}")
    logfile.close()


if __name__ == "__main__":
    main()
