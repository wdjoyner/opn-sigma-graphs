"""Command-line interface: ``python -m sigmagraph N [options]``."""

from __future__ import annotations

import argparse
import json

from .arith import abundancy, classify, factor, sigma
from .graph import SigmaGraph
from .latex import latex_table, tikz


def summary(n: int, internal: bool = False) -> dict:
    G = SigmaGraph(n, extended=not internal)
    out = {
        "n": n,
        "factorization": factor(n),
        "sigma": sigma(n),
        "abundancy": str(abundancy(n)),
        "type": classify(n),
        "vertices": G.vertices,
        "exponents": G.a,
        "boundary": G.boundary(),
        "edges": [[s, r, e] for s, r, e in G.edges()],
        "strongly_connected": G.strongly_connected(),
        "components": G.components(),
        "symmetric": G.is_symmetric(),
    }
    if not internal:
        out["defect"] = G.defect()
        out["defect_identity_holds"] = G.check_defect()
        out["poisson_identity_holds"] = G.check_poisson()
    st = G.arithmetical_structure()
    if st is not None:
        R, S = st
        out["arithmetical_structure"] = {"R": R, "S": S}
        out["critical_group"] = G.critical_group()
    out["edge_types"] = G.classify_edges()
    return out


def main(argv=None) -> None:
    ap = argparse.ArgumentParser(
        prog="sigmagraph",
        description="Valuation-weighted sigma-graph of an integer n > 1.")
    ap.add_argument("n", type=int, help="integer > 1")
    ap.add_argument("--internal", action="store_true",
                    help="use the internal graph (primes dividing n only)")
    ap.add_argument("--json", action="store_true", help="print JSON")
    ap.add_argument("--latex", action="store_true", help="print a LaTeX table")
    ap.add_argument("--tikz", action="store_true", help="print TikZ code")
    args = ap.parse_args(argv)

    G = SigmaGraph(args.n, extended=not args.internal)
    if args.latex:
        print(latex_table(G))
        return
    if args.tikz:
        print(tikz(G))
        return
    info = summary(args.n, args.internal)
    if args.json:
        print(json.dumps(info, indent=2, default=str))
        return
    print(f"n = {args.n} = " + " * ".join(
        f"{p}^{a}" if a > 1 else str(p) for p, a in info["factorization"].items()))
    print(f"sigma(n) = {info['sigma']},  h(n) = {info['abundancy']}  ({info['type']})")
    print(f"vertices: {G.vertices}   boundary: {G.boundary()}")
    print("edges (source -> target : weight):")
    for s, r, e in G.edges():
        print(f"  {s} -> {r} : {e}")
    if not args.internal:
        print(f"defect d(n): {info['defect']}")
    print(f"strongly connected: {info['strongly_connected']}")
    if "arithmetical_structure" in info:
        print(f"arithmetical structure R = {info['arithmetical_structure']['R']}")
        print(f"                       S = {info['arithmetical_structure']['S']}")
        print(f"critical group invariant factors: {info['critical_group'] or 'trivial'}")
