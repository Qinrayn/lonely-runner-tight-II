# Tight instances of the Lonely Runner Conjecture: one-entry modifications and the growth of speeds

Code and data accompanying the paper

> **Tight instances of the Lonely Runner Conjecture: one-entry
> modifications and the growth of speeds**
> Yuhan Zhang (Yangzhou University), 2026.
> Preprint: arXiv:2608.13599.

For a set `V` of `n-1` distinct positive integers put

```
LR(V) = max_t min_{v in V} ||v t||
```

where `||x||` is the distance from `x` to the nearest integer. The Lonely Runner
Conjecture asserts `LR(V) >= 1/n`; `V` is **tight** when `LR(V) = 1/n`. This
repository verifies every computational claim of the paper, in **exact rational
arithmetic**, and reproduces all of its tables and figures.

The paper's main results:
1. **Complete classification** of one-entry modifications `[n-1]_{r -> w}`
   (Theorem 1.2): the proof is purely theoretical, via a key inequality
   `4rI/(2s-I) < r(n-2)` proved for all `r >= 3` using the Jacobsthal function.
2. **Extremal speed** (Theorem 1.3): `max V <= 0.60 n log n + 52 n` for
   every tight one-entry modification, unconditionally; constant `1/2`
   sharp along the primorial sequence `n = p# + 2`.
3. **Deletion values**: a short self-contained proof of Liu's closed
   formula `LR([n-1]\{2}) = 2/(2*floor(n/2)+3)` (Theorem 2.7; the
   general case `r <= (n-1)/2` is Liu's theorem), and
   `LR([n-1]\{r}) > 1/(n-1)` for every `r <= n-2-2*sqrt(n)`
   (Theorem 2.5 and Corollary 2.6).
4. **Two-swap census** (Theorem 6.4): `{1,4,5,6,7,11,13}` (`n = 8`) and the
   Goddyn--Wong doubling member (`n = 74`) are the only tight two-swaps with
   `n <= 140`, `w <= 8n`, over **all** removal pairs; no tight two-swap
   removes a speed `q` with `2 <= q <= n/10` once `n >= max(40, 10q)`
   (Theorem 5.5, covering-density argument, verified by `q_density.py`).
   The doubling family and the
   growth bound under the Goddyn--Wong criterion are verified by
   `k2_gw_284.py` and `theoremB_bound.py`.
5. **Spectrum near the minimum** (Section 7.2): the smallest non-tight
   value of `LR([n-1]_{2 -> w})` is determined for `n = 5..12` over all
   `w >= n` (uniform covering bound, Lemma 7.6).

## Quick start

```bash
pip install -r requirements.txt

# Theorem 1.2: unified proof verification (key inequality, all r >= 3)
python unified_proof.py           # r <= 60 (default, matches paper audit)
python unified_proof.py 200       # extended range

# Section 6: density theorems backing the two-swap exclusions
python r2_density.py              # exact verification of Lemmas 6.1-6.4 + thresholds
python q_density.py               # windows at 1/q, edge reach, inequality thresholds

# Two-swap census: ALL removal pairs, n = 8..140, w <= 8n (~35 min)
python two_swap_census.py --smoke # n = 8..12 positive control first
python two_swap_census.py         # n = 74 k-GW control, then full census
python two_swap_census.py --deep  # deep window: R={2,r'}, n = 9..18, w <= 40n

# Theorem B: extremal speed M(n) ~ (1/2) n log n
python primorial.py               # family table + branch-and-bound certificates
python theoremB_bound.py          # M(n) <= 0.60 n log n + 52 n, n <= 40000
python theta084_check.py          # theta(x) >= 0.84 x for 101 <= x < 523 (exact primorials)
python k2_gw_284.py               # k-GW n=284 certification (60,979 cells)

# Deletion values and the spectrum near the minimum
python verify_liu_unif.py         # Liu's formula (s=1) vs exact LR; Lemma 7.6 thresholds
python spectrum_table2.py         # Table 2: smallest non-tight LR([n-1]_{2->w}), n = 5..12

# Independent cross-checks (legacy searches, exact caps as documented)
python corrected_two_swap_search.py   # union condition, n = 8..24
python two_swap_proof.py              # same-side lemma + global hole, n = 8..17

# Theorem 1.2 proof audit (window/sign lemmas, phi-lemma, threshold chain)
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

### Theorem 1.2 (complete classification)
| file | contents |
|---|---|
| `unified_proof.py` | verifies the key inequality `4rI/(2s-I) < r(n-2)` for all `r >= 3` (Case A/B, Jacobsthal bound). |
| `allcomp_filter.py` | the all-components filter: every component of `U(n,r)` must be covered, not just the longest. |
| `midrange.py` | the mid-range classification as a finite check (`n <= 6r`, `w <= 4rI/(2s-I)`). |
| `verify_stepB.py` | line-by-line numeric audit of the window/sign lemmas. |
| `falsify.py` | self-falsification suite: independent attacks on the weakest joints of the proof. |

### Theorem B (extremal speed), deletion values, and the spectrum
| file | contents |
|---|---|
| `theoremB_bound.py` | `M(n) < 0.60 n log n + 52 n` via Rosser--Schoenfeld; full `s`-scan for `n <= 5000`, `s <= 64` up to `n = 40000`. |
| `theta084_check.py` | the finite branch of the Chebyshev constant: `theta(x) >= 0.84 x` for all integers `101 <= x < 523`, from exact primorials. |
| `primorial.py` | lower bound: the primorial construction attains `w ~ (1/2) n log n`; branch-and-bound certificates. |
| `verify_liu_unif.py` | Liu's closed formula for deletion values (case `s = 1`, `r <= (n-1)/2`) re-derived exactly for `n = 5..24`; the uniform covering bound (Lemma 7.6) thresholds for Table 2. |
| `spectrum_table2.py` | Table 2: the smallest non-tight `LR([n-1]_{2->w})` over `n <= w <= 8n`, exact arithmetic. | |

### Goddyn--Wong multi-acceleration and the two-swap census (Section 6)
| file | contents |
|---|---|
| `two_swap_census.py` | exhaustive census: **all** removal pairs `r1 < r2` (any regime), `n = 8..140`, `w <= 8n`; edge-reach pruning + float screen + exact certification; `--deep` runs the `{2,r'}` window `n = 9..18`, `w <= 40n`. |
| `r2_density.py` | exact verification of the density-theorem chain: `U(n,2)` closed form, isolation lemma, edge reach, and the coverage inequality for odd `n >= 19` / even `n >= 20` (checked for every `19 <= n <= 199`). |
| `k2_gw_284.py` | branch-and-bound certification of the second doubling-family member (`n = 284`, 60,979 cells). |
| `two_swap_23_n8_deep.py` | the `n = 8`, `R = {2,3}` deep window: all 901 candidate pairs `8 <= w1 < w2 <= 320` with `w2 > 64` checked in exact arithmetic; none covers. |
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
| Thm 1.2 (key inequality, all `r >= 3`) | `unified_proof.py` | no violation for `r <= 60` |
| Thm 1.2 (window/sign, no tight mid-range) | `allcomp_filter.py`, `midrange.py`, `verify_stepB.py`, `falsify.py` | `r >= 3`: none tight; `r = 2`: exactly `(5,7),(6,9)` |
| Thm 1.3 (`M(n) < 0.60 n log n + 52 n`) | `theoremB_bound.py`, `theta084_check.py`, `primorial.py` | confirmed for `n <= 40000`; extremal ratio `0.9246` at `n = 212` |
| Thm 2.7 (deletion value for `r = 2`, Liu's formula) | `verify_liu_unif.py` | Liu's theorem (case `s = 1`) re-derived exactly for `n = 5..24`, all `r <= (n-1)/2`: 0 mismatches |
| Lemma 7.6 + Table 2 (spectrum) | `verify_liu_unif.py`, `spectrum_table2.py` | thresholds `12,20,10,14,14,18,18,22` exact; table values certified for all `w >= n` |
| Sec. 6 density theorem (`r=2` exclusion, `n >= 19/20`) | `r2_density.py` | all lemmas verified exactly; inequality holds for every odd `n >= 19`, even `n >= 20` (checked to `n = 199`; crude bound beyond) |
| Doubling family (`n = 74`, `n = 284`) | `two_swap_census.py`, `k2_gw_284.py`, `corrected_two_swap_search.py` | `LR = 1/74` and `LR = 1/284` in exact arithmetic |
| Two-swap census (`n <= 140`, `w <= 8n`, all pairs) | `two_swap_census.py` (+ two legacy cross-checks) | only `{1,4,5,6,7,11,13}` (`n = 8`) and the GW doubling `n = 74`; see `two_swap_census.log` for case counts |
| `{2,r'}` deep window (`n = 9..18`, `w <= 40n`) | `two_swap_census.py --deep` | no tight two-swap |
| `{2,3}` deep window (`n = 8`, `64 < w2 <= 320`) | `two_swap_23_n8_deep.py` | 901 candidate pairs, none covers |
| `k = 3, 4` residual | `complete_finite_verify.py`, `k3_verify_*.py`, `residual_verify.py` | 0 tight beyond the Goddyn--Wong families |

## Environment

Developed with CPython 3.13 on Windows 11; `numpy` is used only for the
pre-screening kernels; `matplotlib` only for the figure. Exact results are
independent of platform and of the presence of `numpy`.

## Licence

Code: MIT (see `LICENSE`).
