# Tight instances of the Lonely Runner Conjecture: complete classification of one-entry modifications, tight multi-swaps, and the growth bound

Code and data accompanying the paper

> **Tight instances of the Lonely Runner Conjecture: complete
> classification of one-entry modifications, tight multi-swaps, and
> the growth bound**
> Yuhan Zhang (Yangzhou University), 2026.

For a set `V` of `n-1` distinct positive integers put

```
LR(V) = max_t min_{v in V} ||v t||
```

where `||x||` is the distance from `x` to the nearest integer. The Lonely Runner
Conjecture asserts `LR(V) >= 1/n`; `V` is **tight** when `LR(V) = 1/n`. This
repository verifies every computational claim of the paper, in **exact rational
arithmetic**, and reproduces all of its tables and figures.

The paper proves four results:
1. **Complete classification** of one-entry modifications `[n-1]_{r -> w}`
   (Theorem A): the proof is purely theoretical, via a key inequality
   `4rI/(2s-I) < r(n-2)` proved for all `r >= 3` using the Jacobsthal function.
2. **Multi-swaps** (Theorems 5.1/5.2): the Goddyn--Wong multi-acceleration
   theorem (their 2006 Thm 3.1) is recalled, and an explicit CRT doubling
   subfamily is given; members `n = 74` (the GW example) and `n = 284`
   are certified in exact arithmetic.
3. **Growth bound** (Theorems B / 6.2): `max V <= 0.60 n log n + 52 n` for
   every tight one-entry modification with `n >= 4` runners, unconditionally;
   constant `1/2` sharp along the primorial sequence `n = p# + 2`.
4. **Two-swap census** (Theorem 7.13): `{1,4,5,6,7,11,13}` (`n = 8`) and the
   Goddyn--Wong doubling instance (`n = 74`) are the only tight two-swaps with `n <= 140`,
   `w <= 8n`, over **all** removal pairs; for removals `{2, r'}` with
   mid-range `r'` the exclusion is theoretical for `n >= 19` (odd) /
   `n >= 20` (even), all `w` (isolation + coverage-density lemmas, verified
   by `r2_density.py`); the same mechanism run at `1/q` excludes two-swaps
   containing **any** removal `2 <= q <= n/10` for `n >= max(40, 10q)`,
   with no restriction on the second removal or the speeds (verified by
   `q_density.py`).
5. **Deletion values** (Theorem 2.8): `LR([n-1] \\ {2}) = 2/(2*floor(n/2)+3)`
   in closed form for all `n >= 5`; `LR([n-1] \\ {r}) >= 1/(n-1)` verified
   for all `n <= 30` (all `r`), with equality iff `r = n-1`.

## Quick start

```bash
pip install -r requirements.txt

# Theorem A: unified proof verification (key inequality, all r >= 3)
python unified_proof.py           # r <= 60 (default, matches paper audit)
python unified_proof.py 200       # extended range

# Section 6: density theorem backing the r=2 two-swap exclusion
python r2_density.py              # exact verification of Lemmas 6.1-6.4 + thresholds

# Section 6.2: general-q density chain (any small removal, any r', all w)
python q_density.py               # windows at 1/q, edge reach, inequality thresholds,
                                  # + deletion-value audits (r=2 closed form, n <= 30)

# Two-swap census: ALL removal pairs, n = 8..140, w <= 8n (~35 min)
python two_swap_census.py --smoke # n = 8..12 positive control first
python two_swap_census.py         # n = 74 GW control, then full census
python two_swap_census.py --deep  # deep window: R={2,r'}, n = 9..18, w <= 40n

# Theorem B: extremal speed M(n) ~ (1/2) n log n
python primorial.py               # family table + branch-and-bound certificates
python theoremB_bound.py          # M(n) <= 0.60 n log n + 52 n, n <= 40000
python k2_gw_284.py               # GW doubling n=284 certification (60,979 cells)

# Independent cross-checks (legacy searches, exact caps as documented)
python corrected_two_swap_search.py   # union condition, n = 8..24
python two_swap_proof.py              # same-side lemma + global hole, n = 8..17

# Theorem A proof audit (window/sign lemmas, phi-lemma, threshold chain)
python allcomp_filter.py
python midrange.py
python verify_stepB.py
python falsify.py

# k = 3, 4 residual verification
python complete_finite_verify.py
python k3_verify_extended.py
python k3_verify_r16_25.py
python residual_verify.py
```

Each script is self-contained and prints a summary of the checks performed.
All exact computations use `fractions.Fraction`; `numpy` is used only for
float pre-screening and is optional for the exact path.

## Correctness design

All reported values of `LR` are **theorems, not estimates**:

1. **Breakpoint certification.** `LR(V)` is attained on the finite set of
   breakpoints `(2a+1)/(2v)`, `b/(v+w)`, `b/(w-v)`. `ml.ml_exact` evaluates
   `f_V` there using `fractions.Fraction` only -- no floating point.
2. **Provably-safe screening.** The large enumerations screen in floats over
   the same breakpoint/interval sets, discarding a candidate only when the
   failure is far outside float noise. In the census (`two_swap_census.py`)
   bad intervals are *expanded* before float subtraction (so every float
   component under-approximates the true uncovered region), and a pair is
   discarded only on a gap exceeding `1e-9`, while true gaps are at least
   `1/(n w1 w2) >= 5.6e-9` at `n = 140`. Every survivor is recomputed
   exactly.
3. **Independent implementations.** The two-swap census uses three algorithms
   (edge-reach census, union filter, remaining-region method) sharing no
   code, with complete agreement; `lrc/verify_independent.py` re-derives
   `LR` along two further algorithmic paths (wall subdivision,
   adaptive branch-and-bound) in exact arithmetic.
4. **Positive controls.** Each search run re-discovers the known tight
   instances (`n = 8` `{11,13}`; `n = 74` `{140,144}`) before scanning.

## Files

### Exact-arithmetic core (from the predecessor repository)
| file | contents |
|---|---|
| `lrc/ml.py` | `ml_exact(V)` -- exact `LR(V)` and a maximiser, via the breakpoint set. |
| `lrc/mu_criterion.py` | `U(n,r)` as an exact union of open rational intervals. |
| `lrc/twoswap_hunt.py` | general two-speed modifications with arbitrary inserted speeds. |
| `lrc/verify_independent.py` | wall subdivision + adaptive branch-and-bound certification. |

### Theorem A (complete classification)
| file | contents |
|---|---|
| `unified_proof.py` | verifies the key inequality `4rI/(2s-I) < r(n-2)` for all `r >= 3` (Case A/B, Jacobsthal bound). |
| `allcomp_filter.py` | the all-components filter: every component of `U(n,r)` must be covered, not just the longest. |
| `midrange.py` | the mid-range classification as a finite check (`n <= 6r`, `w <= 4rI/(2s-I)`). |
| `verify_stepB.py` | line-by-line numeric audit of the window/sign lemmas. |
| `falsify.py` | self-falsification suite: independent attacks on the weakest joints of the proof. |

### Theorem B (extremal speed)
| file | contents |
|---|---|
| `theoremB_bound.py` | `M(n) < 0.60 n log n + 52 n` via Rosser--Schoenfeld; full `s`-scan for `n <= 5000`, `s <= 64` up to `n = 40000`. |
| `primorial.py` | lower bound: the primorial construction attains `w ~ (1/2) n log n`; branch-and-bound certificates. |

### Multi-swap family and two-swap census (Section 5+)
| file | contents |
|---|---|
| `two_swap_census.py` | exhaustive census: **all** removal pairs `r1 < r2` (any regime), `n = 8..140`, `w <= 8n`; edge-reach pruning + float screen + exact certification; `--deep` runs the `{2,r'}` window `n = 9..18`, `w <= 40n`. |
| `r2_density.py` | exact verification of the density-theorem chain: `U(n,2)` closed form, isolation lemma, edge reach, and the coverage inequality for odd `n >= 19` / even `n >= 20` (checked for every `19 <= n <= 199`). |
| `q_density.py` | exact verification of the general-q chain (Section 6.2): windows at `1/q` are exactly the adjacent components of `U(n,{q,r'})` (14,742 instances), generalized edge reach, inequality thresholds `n0(q)`, limit constant `(9-sqrt(57))/12`, plus the deletion-value audits (`r=2` closed form to `n <= 30`, sandwich `LR >= 1/(n-1)` for all `n <= 30`). |
| `k2_gw_284.py` | branch-and-bound certification of the `n = 284` member of the doubling family (60,979 cells). |
| `corrected_two_swap_search.py` | union-condition search, exact caps: mid-range `n = 8..17` `w <= 10n`; `r1=2` `n = 8..16` `w <= 10n`; `{2,3}` `n = 8..24` `w <= 6n`. |
| `two_swap_proof.py` | same-side lemma + global-hole verification, `n = 8..17`, `w <= 8n`. |

### k = 3, 4 residual verification
| file | contents |
|---|---|
| `complete_finite_verify.py` | `k = 3` mid-range removals for `r = 3..15`, `n = 2r+1..4r`. |
| `k3_verify_extended.py` | `k = 3` verification extended to `n = 26..35`. |
| `k3_verify_r16_25.py` | `k = 3` verification for `r = 16..25` (`n` up to `100`). |
| `residual_verify.py` | finite residual: `n < 3r+1` for `r = 3..15`. |

## Main verified statements

| paper result | script | outcome |
|---|---|---|
| Thm A (key inequality, all `r >= 3`) | `unified_proof.py` | no violation for `r <= 60` |
| Thm A (window/sign, no tight mid-range) | `allcomp_filter.py`, `midrange.py`, `verify_stepB.py`, `falsify.py` | `r >= 3`: none tight; `r = 2`: exactly `(5,7),(6,9)` |
| Thm B (`M(n) < 0.60 n log n + 52 n`) | `theoremB_bound.py`, `primorial.py` | confirmed for `n <= 40000`; extremal ratio `0.9246` at `n = 212` |
| Sec. 6 density theorem (`r=2` exclusion, `n >= 19/20`) | `r2_density.py` | all lemmas verified exactly; inequality holds for every odd `n >= 19`, even `n >= 20` (checked to `n = 199`; crude bound beyond) |
| Sec. 6.2 general-q exclusion (`q <= n/10`, any `r'`, all `w`) | `q_density.py` | windows exact (14,742 instances); edge reach 0 violations; inequality holds for `n >= max(40,10q)`; limit constant `0.1208` |
| Thm 2.6 deletion value (`r=2` closed form) | `q_density.py` | `2/(2*floor(n/2)+3)` exact for `5 <= n <= 30`; `LR >= 1/(n-1)` all `(n,r)`, `n <= 30`, equality iff `r = n-1` |
| GW doubling family (`n = 74`, `n = 284`) | `two_swap_census.py`, `k2_gw_284.py`, `corrected_two_swap_search.py` | `LR = 1/74` and `LR = 1/284` in exact arithmetic |
| Two-swap census (`n <= 140`, `w <= 8n`, all pairs) | `two_swap_census.py` (+ two legacy cross-checks) | only `{1,4,5,6,7,11,13}` (`n = 8`) and the GW doubling instance `n = 74`; see `two_swap_census.log` for case counts |
| `{2,r'}` deep window (`n = 9..18`, `w <= 40n`) | `two_swap_census.py --deep` | no tight two-swap |
| `k = 3, 4` residual | `complete_finite_verify.py`, `k3_verify_*.py`, `residual_verify.py` | 0 tight beyond the GW families |

## Environment

Developed with CPython 3.13 on Windows 11; `numpy` is used only for the
pre-screening kernels; `matplotlib` only for the figure. Exact results are
independent of platform and of the presence of `numpy`.

## Licence

Code: MIT (see `LICENSE`).
