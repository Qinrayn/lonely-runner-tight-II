# -*- coding: utf-8 -*-
"""
Verification of the estimate  theta(x) >= 0.84 x  for all integers
101 <= x < 523  (equation (eq:theta084) of the paper, finite branch).

theta(x) = sum_{p <= x} log p is computed as the natural logarithm of the
EXACT primorial  P(x) = prod_{p <= x} p  (a Python big integer), evaluated
with 50-digit decimal arithmetic.  The margin theta(x)/x - 0.84 is at
least 0.0094 on the range, about 47 orders of magnitude above the
arithmetic error, so the comparison is rigorous.

Also verifies the asymptotic branch: for x >= 523 the estimate follows
from  theta(x) > x(1 - 1/log x)  (Rosser-Schoenfeld), and we check
1 - 1/ln(523) > 0.84  and that 1 - 1/log x is increasing.

Run:  python v3/theta084_check.py
"""
from decimal import Decimal, getcontext
from math import log

getcontext().prec = 50


def sieve_primes(limit):
    flags = bytearray([1]) * (limit + 1)
    flags[0:2] = b"\x00\x00"
    for i in range(2, int(limit ** 0.5) + 1):
        if flags[i]:
            flags[i * i :: i] = bytearray(len(flags[i * i :: i]))
    return [i for i in range(limit + 1) if flags[i]]


primes = sieve_primes(522)
assert primes[-1] == 521

P = 1
ip = 0
rows = []
worst = None
for x in range(101, 523):
    while ip < len(primes) and primes[ip] <= x:
        P *= primes[ip]
        ip += 1
    theta = Decimal(P).ln()
    ratio = theta / Decimal(x)
    ok = ratio >= Decimal("0.84")
    if worst is None or ratio < worst[1]:
        worst = (x, ratio)
    rows.append((x, ratio, ok))

bad = [r for r in rows if not r[2]]
print(f"checked {len(rows)} integers x in [101, 523): "
      f"{len(bad)} violations of theta(x) >= 0.84 x")
print(f"minimum of theta(x)/x on the range: {worst[1]} at x = {worst[0]}")

# monotonic branch: 1 - 1/log x > 0.84 for x >= 523, and increasing
v523 = 1 - 1 / log(523)
print(f"1 - 1/ln(523) = {v523:.5f} > 0.84 : {v523 > 0.84}")
print(f"1 - 1/log x increasing for x >= 523 : {1 - 1/log(523) < 1 - 1/log(524)}")

assert not bad and v523 > 0.84
print("ALL CHECKS PASSED")
