# -*- coding: utf-8 -*-
# Reproduce the n=8, R={2,3} deep-window exclusion (64 < max(w) <= 40n=320)
# claimed in the proof of the R={2,3} two-swap classification theorem.
import math
from fractions import Fraction as Fr

n = 8
S = [1, 4, 5, 6, 7]        # retained speeds (removed {2,3})
th = Fr(1, 8)

# exact breakpoints of U(8,{2,3}): ||s t|| = 1/8  <=>  t = (8j +/- 1)/(8s)
bps = {Fr(0), Fr(1)}
for s in S:
    for j in range(0, s + 1):
        for num in (8*j - 1, 8*j + 1):
            t = Fr(num, 8*s)
            if 0 <= t <= 1:
                bps.add(t)
pts = sorted(bps)

def inU(t):
    for s in S:
        x = (s*t) % 1
        if min(x, 1-x) <= th:
            return False
    return True

comps = [(a, b) for a, b in zip(pts, pts[1:]) if b > a and inU(Fr(a+b, 2))]
print("U(8,{2,3}): %d components" % len(comps))
for a, b in comps:
    print("   [%.7f, %.7f]  length %.7f" % (float(a), float(b), float(b-a)))

def covers(w1, w2):
    # closed-interval coverage of every closed component by bad sets of w1,w2
    for (a, b) in comps:
        ivs = []
        for w in (w1, w2):
            for j in range(math.ceil(a*w - th), math.floor(b*w + th) + 1):
                c = Fr(j, w); h = Fr(1, 8*w)
                if c + h >= a and c - h <= b:
                    ivs.append((c - h, c + h))
        ivs.sort()
        cur = a; done = False
        for L, R in ivs:
            if L > cur:
                return False
            if R > cur:
                cur = R
            if cur >= b:
                done = True; break
        if not done:
            return False
    return True

print("validation: {11,13} covers U (expect True):", covers(11, 13))
print("validation: {11,15} covers:", covers(11, 15))

def act(w, a):   # some bad interval of w contains the component end a
    return math.ceil(a*w - th) <= math.floor(a*w + th)

WLO, WHI = 8, 320
K = len(comps)
act_tab = {w: [act(w, a) for (a, b) in comps] for w in range(WLO, WHI + 1)}

def mask_pass(w1, w2):
    a1, a2 = act_tab[w1], act_tab[w2]
    return all(a1[k] or a2[k] for k in range(K))

# domain claimed in the proof: pairs with at least one speed in (64, 320]
cnt = 0; covering = []
for w2 in range(65, WHI + 1):
    for w1 in range(WLO, w2):
        if mask_pass(w1, w2):
            cnt += 1
            if covers(w1, w2):
                covering.append((w1, w2))
print("candidate pairs (8<=w1<w2<=320, w2>64) after edge-reach mask:", cnt)
print("covering pairs among them:", covering)

# reference: full range up to 320 (incl. census-covered both<=64)
cnt2 = 0; cov2 = []
for w2 in range(WLO+1, WHI + 1):
    for w1 in range(WLO, w2):
        if mask_pass(w1, w2):
            cnt2 += 1
            if covers(w1, w2):
                cov2.append((w1, w2))
print("reference: candidates over ALL pairs <=320:", cnt2, "; covering:", cov2)

# if any covering pair with w2>64 exists, compute exact LR to test tightness
def LRval(V):
    V = sorted(set(V)); cands = {Fr(0)}
    for v in V:
        for k in range(v): cands.add(Fr(2*k+1, 2*v))
    for i in range(len(V)):
        for j in range(i+1, len(V)):
            u, w = V[i], V[j]
            for k in range(1, w-u): cands.add(Fr(k, w-u))
            for k in range(1, u+w): cands.add(Fr(k, u+w))
    def nrm(v, t):
        x = (v*t) % 1
        return min(x, 1-x)
    return max(min(nrm(v, t) for v in V) for t in cands)

for (w1, w2) in covering[:5]:
    val = LRval([1, 4, 5, 6, 7, w1, w2])
    print("LR({1,4,5,6,7,%d,%d}) = %s  (1/8 = %s)  tight: %s" % (w1, w2, val, Fr(1,8), val == Fr(1,8)))

print("new enumeration total would be: 820404 + %d = %d" % (cnt, 820404 + cnt))
