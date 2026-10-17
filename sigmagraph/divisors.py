"""Divisor graphs and Hasse graphs (Section 2.7 of the notes).

* The divisor lattice Div(n) is the poset of divisors under divisibility.
* The divisor graph D(n) is its comparability graph.
* The Hasse graph H(n) is the undirected graph of covering relations
  d -- d*p, isomorphic to a Cartesian product of paths.
"""

from __future__ import annotations

from fractions import Fraction
from math import prod

from sympy import divisors, nextprime

from .arith import factor, sigma, abundancy


def divisor_graph_edges(n: int) -> list[tuple[int, int]]:
    """Edges {d, d'} of D(n) with d | d', d < d'."""
    ds = divisors(n)
    return [(d, e) for i, d in enumerate(ds) for e in ds[i + 1:] if e % d == 0]


def hasse_graph_edges(n: int) -> list[tuple[int, int, int]]:
    """Edges (d, d*p, p) of H(n), colored by the prime p."""
    ps = sorted(factor(n))
    return [(d, d * p, p) for d in divisors(n) for p in ps if n % (d * p) == 0]


def exponent_pattern(n: int) -> tuple[int, ...]:
    """Sorted exponents; D(n) and H(n) depend only on this (unlabeled)."""
    return tuple(sorted(factor(n).values()))


def partition_function(n: int) -> int:
    """Z(n) = sum of the vertex labels of the labeled Hasse graph = sigma(n)."""
    return sum(divisors(n))


def edge_weighted_abundancy(n: int) -> Fraction:
    """h(n) computed from edge weights 1/p on H(n).

    The weight of an upward path from d to n is the product of 1/p over its
    edges, which is d/n for every path.  Summing over vertices gives h(n).
    The path used here climbs one prime at a time in increasing order.
    """
    fn = factor(n)
    total = Fraction(0)
    for d in divisors(n):
        fd = factor(d) if d > 1 else {}
        w = Fraction(1)
        for p, a in sorted(fn.items()):
            w *= Fraction(1, p) ** (a - fd.get(p, 0))
        total += w
    return total


def top_weighted_degree(n: int) -> Fraction:
    """Weighted degree of the top vertex n in D(n) with weights d/d'.

    Equals h(n) - 1, so n is perfect iff this is 1.
    """
    return sum((Fraction(d, n) for d in divisors(n) if d < n), Fraction(0))


def deficient_twin(n: int) -> int:
    """An odd deficient number with the same exponent pattern as n.

    Replaces the k primes of n by k distinct primes > 2k, as in the
    proposition that unlabeled divisor graphs cannot detect perfectness.
    """
    exps = sorted(factor(n).values())
    k = len(exps)
    q, primes = 2 * k, []
    while len(primes) < k:
        q = nextprime(q)
        primes.append(int(q))
    m = prod(p ** a for p, a in zip(primes, exps))
    assert abundancy(m) < 2
    return m
