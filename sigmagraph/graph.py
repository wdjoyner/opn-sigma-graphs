"""The sigma-graph of an integer (Section 2 of the notes).

For n > 1 the *extended* sigma-graph has vertex set the primes dividing
n * sigma(n); the weight of s -> r is v_r(sigma(s**a_s)) where a_s = v_s(n)
(a_s = 0 for "boundary" primes that divide sigma(n) but not n).  The
*internal* sigma-graph is the induced subgraph on the primes dividing n.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from fractions import Fraction

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

from .arith import abundancy, edge_type, factor, local_sigma, sigma, valuation


@dataclass
class SigmaGraph:
    """Weighted digraph attached to an integer n > 1.

    Attributes
    ----------
    n         : the integer
    vertices  : sorted list of primes
    a         : exponent v_s(n) of each vertex (0 on boundary vertices)
    W         : dict (s, r) -> positive weight
    """

    n: int
    extended: bool = True
    vertices: list[int] = field(init=False)
    a: dict[int, int] = field(init=False)
    W: dict[tuple[int, int], int] = field(init=False)

    def __post_init__(self):
        if self.n < 2:
            raise ValueError("n must be > 1")
        fn = factor(self.n)
        primes = set(fn)
        if self.extended:
            primes |= set(factor(sigma(self.n)))
        self.vertices = sorted(primes)
        self.a = {p: fn.get(p, 0) for p in self.vertices}
        self.W = {}
        for s in self.vertices:
            if self.a[s] == 0:
                continue
            for r, e in factor(local_sigma(s, self.a[s])).items():
                if r in self.a:
                    self.W[(s, r)] = e

    # ---- basic data -------------------------------------------------
    def w(self, s: int, r: int) -> int:
        return self.W.get((s, r), 0)

    def edges(self) -> list[tuple[int, int, int]]:
        return sorted((s, r, e) for (s, r), e in self.W.items())

    def boundary(self) -> list[int]:
        """Vertices with exponent 0 (no outgoing edges)."""
        return [p for p in self.vertices if self.a[p] == 0]

    def in_degree(self, r: int) -> int:
        return sum(self.w(s, r) for s in self.vertices)

    def out_degree(self, s: int) -> int:
        return sum(self.w(s, r) for r in self.vertices)

    def matrix(self) -> Matrix:
        """Weight matrix W (rows = sources, columns = targets)."""
        V = self.vertices
        return Matrix(len(V), len(V), lambda i, j: self.w(V[i], V[j]))

    # ---- defect and Poisson equation ---------------------------------
    def defect(self) -> dict[int, int]:
        """d_r(n) = in-degree(r) - a_r.  (Proposition 'Valuation defect'.)"""
        return {r: self.in_degree(r) - self.a[r] for r in self.vertices}

    def check_defect(self) -> bool:
        """Verify h(n) = prod r**d_r exactly (extended graph only)."""
        if not self.extended:
            raise ValueError("defect identity needs the extended graph")
        h = Fraction(1)
        for r, d in self.defect().items():
            h *= Fraction(r) ** d
        return h == abundancy(self.n)

    def poisson_operator(self) -> Matrix:
        """L = diag(a) - W."""
        V = self.vertices
        return Matrix.diag(*[self.a[v] for v in V]) - self.matrix()

    def eta(self) -> list[float]:
        """Local log-abundancies eta_s = log h(s**a_s)."""
        return [math.log(local_sigma(s, self.a[s])) - self.a[s] * math.log(s)
                for s in self.vertices]

    def check_poisson(self, tol: float = 1e-9) -> bool:
        """Check L*ell = -eta and 1^T L = -d(n)^T (Proposition 'Poisson')."""
        V = self.vertices
        L = self.poisson_operator()
        ell = [math.log(v) for v in V]
        row_ok = all(abs(sum(float(L[i, j]) * ell[j] for j in range(len(V)))
                         + e) < tol for i, e in enumerate(self.eta()))
        d = self.defect()
        col_ok = all(sum(L[i, j] for i in range(len(V))) == -d[V[j]]
                     for j in range(len(V)))
        return row_ok and col_ok

    # ---- connectivity ------------------------------------------------
    def _reach(self, start: int, forward: bool = True) -> set[int]:
        seen, stack = {start}, [start]
        while stack:
            x = stack.pop()
            for y in self.vertices:
                e = self.w(x, y) if forward else self.w(y, x)
                if e and y not in seen:
                    seen.add(y)
                    stack.append(y)
        return seen

    def strongly_connected(self) -> bool:
        if len(self.vertices) <= 1:
            return True
        v0 = self.vertices[0]
        full = set(self.vertices)
        return self._reach(v0, True) == full == self._reach(v0, False)

    def components(self) -> list[list[int]]:
        """Strongly connected components."""
        left, comps = set(self.vertices), []
        while left:
            v = min(left)
            comp = self._reach(v, True) & self._reach(v, False)
            comps.append(sorted(comp))
            left -= comp
        return comps

    def is_symmetric(self) -> bool:
        return all(self.w(s, r) == self.w(r, s)
                   for s in self.vertices for r in self.vertices)

    # ---- edge classification -----------------------------------------
    def classify_edges(self) -> list[dict]:
        """Order/LTE type of every edge between odd primes."""
        out = []
        for s, r, e in self.edges():
            if s == 2 or r == 2:
                out.append({"source": s, "target": r, "weight": e,
                            "case": "involves 2"})
                continue
            info = edge_type(s, self.a[s], r)
            info.update(source=s, target=r)
            out.append(info)
        return out

    # ---- arithmetical structure --------------------------------------
    def arithmetical_structure(self):
        """Canonical out-structure (R, S) with S = in-degrees.

        Solves (diag(S) - W) R = 0 for a primitive positive integer R.
        Returns None if the graph is not strongly connected or some
        in-degree is 0.  For a k-perfect n, S = a + v(k).
        """
        V = self.vertices
        S = [self.in_degree(v) for v in V]
        if min(S) == 0 or not self.strongly_connected():
            return None
        M = Matrix.diag(*S) - self.matrix()
        ns = M.nullspace()
        if len(ns) != 1:
            return None
        vec = ns[0]
        den = math.lcm(*[int(Fraction(str(x)).denominator) for x in vec])
        ints = [int(Fraction(str(x)) * den) for x in vec]
        g = math.gcd(*ints)
        R = [x // g for x in ints]
        if R[0] < 0:
            R = [-x for x in R]
        if min(R) <= 0:
            return None
        return dict(zip(V, R)), dict(zip(V, S))

    def critical_group(self) -> list[int]:
        """Invariant factors > 1 of the torsion of coker(diag(S) - W)."""
        V = self.vertices
        S = [self.in_degree(v) for v in V]
        M = Matrix.diag(*S) - self.matrix()
        snf = smith_normal_form(M, domain=ZZ)
        diag = [abs(int(snf[i, i])) for i in range(len(V))]
        return [d for d in diag if d > 1]

    # ---- quadratic-form identity --------------------------------------
    def energy_identity(self) -> tuple[float, float]:
        """Both sides of the quadratic-form identity of the notes:

            sum_{s,r} w(s,r) log^2(s/r) + sum_s (2 a_s - out_s - in_s) log^2 s
              = -2 sum_s eta_s log s.

        Since in_s = a_s + d_s, the coefficient equals a_s - out_s - d_s.
        It holds for every n on the extended graph.
        """
        V = self.vertices
        lg = {v: math.log(v) for v in V}
        lhs = sum(e * (lg[s] - lg[r]) ** 2 for (s, r), e in self.W.items())
        lhs += sum((2 * self.a[s] - self.out_degree(s) - self.in_degree(s))
                   * lg[s] ** 2 for s in V)
        rhs = -2 * sum(e * lg[s] for e, s in zip(self.eta(), V))
        return lhs, rhs
