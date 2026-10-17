"""Regression tests: every numerical example in the notes.

Section and example names refer to notes/AI_OPN-conj_sigma-graphs.tex.
"""

from fractions import Fraction

import pytest

from sigmagraph import (DESCARTES_SPOOF, MULTIPERFECT, SigmaGraph, abundancy,
                        classify, cyclotomic_layers, deficient_twin,
                        divisor_graph_edges, edge_type,
                        edge_weighted_abundancy, exponent_pattern, factor,
                        hasse_graph_edges, local_sigma, order,
                        partition_function, primitive_prime_divisors, sigma,
                        sigma0, sigma_r, spoof_sigma, top_weighted_degree,
                        valuation)


# ---------------------------------------------------------------- Section 1
def test_divisor_functions_of_12():
    assert sigma0(12) == 6 and sigma(12) == 28
    assert sigma_r(12, 0) == 6 and sigma_r(12, 1) == 28
    assert sigma_r(12, 2) == 210
    assert sigma_r(12, -1) == Fraction(7, 3)


@pytest.mark.parametrize("n", [12, 28, 45, 360, 945])
def test_reflection(n):
    for r in (1, 2, 3):
        assert sigma_r(n, -r) == Fraction(1, n ** r) * sigma_r(n, r)


@pytest.mark.parametrize("n,h,kind", [
    (6, Fraction(2), "perfect"), (28, Fraction(2), "perfect"),
    (12, Fraction(7, 3), "abundant"), (15, Fraction(8, 5), "deficient"),
    (45, Fraction(26, 15), "deficient"), (945, Fraction(128, 63), "abundant"),
    (120, Fraction(3), "3-perfect"),
])
def test_catalog(n, h, kind):
    assert abundancy(n) == h and classify(n) == kind


def test_euclid_euler():
    for p in (2, 3, 5, 7, 13):
        N = 2 ** (p - 1) * (2 ** p - 1)
        assert sigma(N) == 2 * N


def test_descartes_spoof():
    D = 3**2 * 7**2 * 11**2 * 13**2 * 22021
    assert D == 198585576189
    assert factor(22021) == {19: 2, 61: 1}
    assert spoof_sigma(DESCARTES_SPOOF) == 2 * D      # formal identity
    assert sigma(D) != 2 * D                           # not genuinely perfect
    assert [local_sigma(p, 2) for p in (3, 7, 11, 13)] == [13, 57, 133, 183]
    assert 22021 + 1 == 2 * 7 * 11**2 * 13


def test_orders():
    assert order(13, 3) == 1
    assert order(13, 61) == 3 and order(61, 13) == 3
    assert order(3, 11) == 5


def test_cyclotomic_layers():
    assert cyclotomic_layers(3, 8) == {3: 13, 9: 757}
    assert 13 * 757 == local_sigma(3, 8) == 9841
    assert cyclotomic_layers(3, 12) == {13: 797161}
    assert factor(797161) == {797161: 1}
    assert factor(cyclotomic_layers(67, 2)[3]) == {3: 1, 7: 2, 31: 1}


def test_primitive_divisor_13_cubed():
    assert factor(13**3 - 1) == {2: 2, 3: 2, 61: 1}
    assert primitive_prime_divisors(13, 3) == [61]
    assert 61 % 6 == 1


# order/LTE classification: the three sources and the Rodgers counterexamples
def test_order_one_edge_13_to_3():
    assert factor(local_sigma(13, 2)) == {3: 1, 61: 1}
    t = edge_type(13, 2, 3)
    assert t["case"] == "order-one" and t["weight"] == t["predicted"] == 1


def test_wieferich_type_3_to_11():
    assert local_sigma(3, 4) == 121
    t = edge_type(3, 4, 11)
    assert t["case"] == "order-d" and t["base"] == 2 and t["wieferich"]
    assert t["weight"] == t["predicted"] == 2


def test_multiplier_3_to_13():
    t = edge_type(3, 38, 13)
    assert t["case"] == "order-d" and t["base"] == 1 and t["multiplier"] == 1
    assert t["weight"] == t["predicted"] == 2 == valuation(13, local_sigma(3, 38))


def test_order_one_weight_two_7_to_3():
    assert factor(local_sigma(7, 8)) == {3: 2, 19: 1, 37: 1, 1063: 1}
    t = edge_type(7, 8, 3)
    assert t["case"] == "order-one" and t["weight"] == 2
    # the (wrong) two-term formula with d = 1 would give 3
    assert valuation(3, 7 - 1) + valuation(3, 9) == 3


def test_order_one_edges_are_common():
    from sympy import primerange
    bad = tot = 0
    for q in primerange(3, 200):
        for b in (1, 2, 3, 4):
            for r in factor(local_sigma(q, 2 * b)):
                tot += 1
                bad += order(q, r) == 1
    assert (bad, tot) == (58, 497)


def test_primitive_two_cycle_13_61():
    assert factor(local_sigma(61, 2)) == {3: 1, 13: 1, 97: 1}
    assert 61 % 13 == 9 and pow(9, 3, 13) == 1


@pytest.mark.parametrize("q", [7, 13, 19, 31, 37, 43])
def test_q_one_mod_three_gives_three(q):
    assert local_sigma(q, 2) % 3 == 0 and local_sigma(q, 2) % 9 == 3


# ---------------------------------------------------------------- Section 2
def test_sigma_graph_6_and_28():
    G = SigmaGraph(6)
    assert G.edges() == [(2, 3, 1), (3, 2, 2)]
    assert G.boundary() == []
    G = SigmaGraph(28)
    assert G.edges() == [(2, 7, 1), (7, 2, 3)]


@pytest.mark.parametrize("p", [2, 3, 5, 7, 13])
def test_even_perfect_two_cycle(p):
    M = 2**p - 1
    G = SigmaGraph(2 ** (p - 1) * M)
    assert G.edges() == [(2, M, 1), (M, 2, p)]
    R, S = G.arithmetical_structure()
    assert R == {2: 1, M: p} and S == {2: p, M: 1}
    assert G.critical_group() == []


def test_sigma_graph_15():
    G = SigmaGraph(15)
    assert G.edges() == [(3, 2, 2), (5, 2, 1), (5, 3, 1)]
    assert G.boundary() == [2]
    assert not SigmaGraph(15, extended=False).strongly_connected()


def test_sigma_graph_12():
    G = SigmaGraph(12)
    assert G.edges() == [(2, 7, 1), (3, 2, 2)]
    assert G.boundary() == [7]


def test_sigma_graph_945():
    G = SigmaGraph(945)
    assert G.edges() == [(3, 2, 3), (3, 5, 1), (5, 2, 1), (5, 3, 1), (7, 2, 3)]


def test_sigma_graph_120():
    G = SigmaGraph(120)
    assert G.edges() == [(2, 3, 1), (2, 5, 1), (3, 2, 2), (5, 2, 1), (5, 3, 1)]
    assert [G.in_degree(v) for v in G.vertices] == [3, 2, 1]


@pytest.mark.parametrize("n,b,a,d", [
    (6, (2, 1), (1, 1), (1, 0)),
    (28, (3, 1), (2, 1), (1, 0)),
    (15, (3, 1, 0), (0, 1, 1), (3, 0, -1)),
    (12, (2, 0, 1), (2, 1, 0), (0, -1, 1)),
    (945, (7, 1, 1, 0), (0, 3, 1, 1), (7, -2, 0, -1)),
    (120, (3, 2, 1), (3, 1, 1), (0, 1, 0)),
])
def test_defect_table(n, b, a, d):
    G = SigmaGraph(n)
    V = G.vertices
    assert tuple(G.in_degree(v) for v in V) == b
    assert tuple(G.a[v] for v in V) == a
    assert tuple(G.defect()[v] for v in V) == d
    assert G.check_defect()


@pytest.mark.parametrize("n", [6, 12, 15, 28, 45, 120, 945, 360, 2178540])
def test_poisson_identity(n):
    assert SigmaGraph(n).check_poisson()


def test_poisson_operator_120():
    L = SigmaGraph(120).poisson_operator()
    assert L.tolist() == [[3, -1, -1], [-2, 1, 0], [-1, -1, 1]]


def test_poisson_operator_6():
    assert SigmaGraph(6).poisson_operator().tolist() == [[1, -1], [-2, 1]]


@pytest.mark.parametrize("n,R,S,K", [
    (120, (1, 1, 2), (3, 2, 1), [2]),
    (672, (1, 1, 3), (5, 2, 1), [2]),
    (30240, (1, 2, 3, 3), (7, 3, 1, 1), [2]),
    (32760, (1, 2, 3, 3, 4), (5, 2, 1, 1, 1), [2]),
    (2178540, (2, 5, 7, 8, 10, 11), (4, 2, 1, 2, 1, 1), []),
])
def test_arithmetical_structures(n, R, S, K):
    G = SigmaGraph(n)
    Rg, Sg = G.arithmetical_structure()
    V = G.vertices
    assert tuple(Rg[v] for v in V) == R
    assert tuple(Sg[v] for v in V) == S
    assert G.critical_group() == K


def test_factorization_2178540():
    assert factor(2178540) == {2: 2, 3: 2, 5: 1, 7: 2, 13: 1, 19: 1}


@pytest.mark.parametrize("n", MULTIPERFECT)
def test_multiperfect_list(n):
    assert sigma(n) % n == 0


@pytest.mark.parametrize("n", MULTIPERFECT[:14])
def test_multiperfect_properties(n):
    G = SigmaGraph(n)
    k = sigma(n) // n
    assert G.check_defect()
    assert all(G.defect()[r] == valuation(r, k) for r in G.vertices)
    assert not G.is_symmetric()                    # max-principle theorem
    lhs, rhs = G.energy_identity()
    assert abs(lhs - rhs) < 1e-9


def test_energy_identity_120_values():
    lhs, rhs = SigmaGraph(120).energy_identity()
    assert abs(lhs + 2.0904092) < 1e-6 and abs(rhs + 2.0904092) < 1e-6


@pytest.mark.parametrize("n", [6, 28, 496, 8128])
def test_perfect_graphs_strongly_connected(n):
    assert SigmaGraph(n).strongly_connected()


# symmetric (sigma-symmetric) examples
def test_symmetric_pairs():
    for n, (s, r) in [(3**2 * 13**2, (3, 13)), (13**2 * 61**2, (13, 61))]:
        G = SigmaGraph(n, extended=False)
        assert G.w(s, r) == G.w(r, s) == 1


def test_symmetric_closed_135():
    G = SigmaGraph(135)        # 3^3 * 5
    assert sigma(135) == 2**4 * 3 * 5
    odd = [v for v in G.vertices if v != 2]
    assert odd == [3, 5] and G.w(3, 5) == G.w(5, 3) == 1
    assert abundancy(135) == Fraction(16, 9)


# divisor graphs and Hasse graphs
def test_divisor_graph_6():
    assert len(divisor_graph_edges(6)) == 5        # K4 minus the edge {2,3}
    assert (2, 3) not in divisor_graph_edges(6)


def test_hasse_graph_28_and_945():
    assert len(hasse_graph_edges(28)) == 7          # 2x3 grid
    # P4 x P2 x P2: 3 edges in each of 4 copies of P4, plus 8 + 8 others
    assert len(hasse_graph_edges(945)) == 3 * 4 + 8 + 8
    assert sigma0(945) == 16
    assert partition_function(945) == sigma(945) == 1920 == 40 * 6 * 8


def test_one_graph_three_kinds():
    assert exponent_pattern(12) == exponent_pattern(28) == exponent_pattern(45)
    assert exponent_pattern(6) == exponent_pattern(15) == exponent_pattern(35)


def test_edge_weighted_abundancy():
    assert edge_weighted_abundancy(12) == Fraction(7, 3)
    for n in (6, 28, 945, 120):
        assert edge_weighted_abundancy(n) == abundancy(n)
        assert top_weighted_degree(n) == abundancy(n) - 1


@pytest.mark.parametrize("n", [945, 3**4 * 5**2 * 7, 105, 3**2 * 5 * 7 * 11])
def test_deficient_twin(n):
    m = deficient_twin(n)
    assert m % 2 == 1 and abundancy(m) < 2
    assert exponent_pattern(m) == exponent_pattern(n)
