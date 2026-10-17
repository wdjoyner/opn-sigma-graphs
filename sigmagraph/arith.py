"""Elementary number theory used by the sigma-graph toolkit.

Everything here is exact integer or rational arithmetic.  The functions
mirror Section 1 of the notes: divisor functions, abundancy, valuations,
multiplicative orders, and the order/LTE classification of the prime
factors of a local divisor sum.
"""

from __future__ import annotations

from fractions import Fraction
from functools import reduce
from math import gcd, prod

from sympy import factorint, divisors, n_order, cyclotomic_poly
from sympy.abc import x as _x


def factor(n: int) -> dict[int, int]:
    """Prime factorization of n > 0 as {prime: exponent}."""
    if n < 1:
        raise ValueError("n must be a positive integer")
    return {int(p): int(e) for p, e in factorint(n).items()}


def sigma_r(n: int, r: int = 1) -> Fraction:
    """Divisor function sigma_r(n) = sum of d**r over divisors d of n.

    Returns a Fraction so that negative r is exact.
    """
    return sum((Fraction(d) ** r for d in divisors(n)), Fraction(0))


def sigma(n: int) -> int:
    """Sum of divisors sigma(n) = sigma_1(n), computed multiplicatively."""
    return prod((p ** (a + 1) - 1) // (p - 1) for p, a in factor(n).items())


def sigma0(n: int) -> int:
    """Number of divisors sigma_0(n)."""
    return prod(a + 1 for a in factor(n).values())


def local_sigma(p: int, a: int) -> int:
    """Local divisor sum sigma(p**a) = 1 + p + ... + p**a."""
    return (p ** (a + 1) - 1) // (p - 1) if a > 0 else 1


def abundancy(n: int) -> Fraction:
    """Abundancy h(n) = sigma(n)/n."""
    return Fraction(sigma(n), n)


def classify(n: int) -> str:
    """'perfect', 'abundant', 'deficient', or 'k-perfect' (k >= 3)."""
    h = abundancy(n)
    if h.denominator == 1 and h >= 3:
        return f"{h.numerator}-perfect"
    if h == 2:
        return "perfect"
    return "abundant" if h > 2 else "deficient"


def valuation(r: int, x: int | Fraction) -> int:
    """r-adic valuation of a nonzero integer or rational."""
    x = Fraction(x)
    if x == 0:
        raise ValueError("valuation of 0 is undefined")
    def v(m: int) -> int:
        m, k = abs(m), 0
        while m % r == 0:
            m //= r
            k += 1
        return k
    return v(x.numerator) - v(x.denominator)


def omega(n: int) -> int:
    """Number of distinct prime factors."""
    return len(factor(n))


def big_omega(n: int) -> int:
    """Number of prime factors counted with multiplicity."""
    return sum(factor(n).values())


def order(s: int, r: int) -> int:
    """Multiplicative order of s modulo the prime r (requires r not | s)."""
    return int(n_order(s, r))


def cyclotomic_value(d: int, s: int) -> int:
    """Phi_d(s)."""
    return int(cyclotomic_poly(d, _x).subs(_x, s))


def cyclotomic_layers(s: int, a: int) -> dict[int, int]:
    """The layers Phi_d(s), d | a+1, d > 1, whose product is sigma(s**a)."""
    return {d: cyclotomic_value(d, s) for d in divisors(a + 1) if d > 1}


def edge_type(s: int, a: int, r: int) -> dict:
    """Order/LTE classification of v_r(sigma(s**a)) for odd primes s != r.

    Returns a dict with keys
      weight      : v_r(sigma(s**a))
      order       : d = ord_r(s)
      case        : 'none' | 'order-one' | 'order-d'
      base        : v_r(s**d - 1)          (order-d case only)
      multiplier  : v_r((a+1)/d)           (order-d case only)
      wieferich   : base >= 2              (order-d case only)
    matching Proposition 'Order/LTE classification' of the notes.
    """
    if s == r or s % 2 == 0 or r % 2 == 0:
        raise ValueError("s and r must be distinct odd primes")
    n = a + 1
    d = order(s, r)
    out = {"order": d, "weight": valuation(r, local_sigma(s, a))}
    if n % d != 0:
        out["case"] = "none"
    elif d == 1:
        out["case"] = "order-one"
        out["predicted"] = valuation(r, n)
    else:
        base = valuation(r, s ** d - 1)
        mult = valuation(r, n // d)
        out.update(case="order-d", base=base, multiplier=mult,
                   wieferich=base >= 2, predicted=base + mult)
    return out


def primitive_prime_divisors(s: int, d: int) -> list[int]:
    """Primes l with ord_l(s) = d (the primitive prime divisors of s**d - 1)."""
    return sorted(p for p in factor(s ** d - 1)
                  if s % p != 0 and order(s, p) == d)


def spoof_sigma(factorization: dict[int, int]) -> int:
    """Formal divisor sum, treating every key as if it were prime.

    For a genuine factorization this equals sigma(n).  For a spoof such as
    Descartes' {3:2, 7:2, 11:2, 13:2, 22021:1} it does not.
    """
    return prod(local_sigma(p, a) for p, a in factorization.items())


def lcm_list(xs) -> int:
    return reduce(lambda u, v: u * v // gcd(u, v), xs, 1)
