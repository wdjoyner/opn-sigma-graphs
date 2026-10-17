"""Survey of the sigma-graphs of small multiperfect numbers.

Usage:  python examples/survey.py [--csv out.csv] [--count 14]

For each n it prints k = sigma(n)/n, the number of vertices, whether the
extended sigma-graph is strongly connected, the canonical arithmetical
structure R (with S = in-degrees), and the critical group.
"""

import argparse
import csv
import sys

from sigmagraph import MULTIPERFECT, SigmaGraph, sigma


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", help="write results to this CSV file")
    ap.add_argument("--count", type=int, default=14,
                    help="how many multiperfect numbers to survey")
    args = ap.parse_args()

    rows = []
    for n in MULTIPERFECT[: args.count]:
        G = SigmaGraph(n)
        st = G.arithmetical_structure()
        R = [st[0][v] for v in G.vertices] if st else None
        rows.append({
            "n": n,
            "k": sigma(n) // n,
            "vertices": G.vertices,
            "strongly_connected": G.strongly_connected(),
            "R": R,
            "S": [G.in_degree(v) for v in G.vertices],
            "critical_group": G.critical_group() if st else None,
        })

    for r in rows:
        print(f"{r['n']:>14}  k={r['k']}  V={r['vertices']}  "
              f"SC={r['strongly_connected']}  R={r['R']}  K={r['critical_group']}")

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        print(f"wrote {args.csv}", file=sys.stderr)


if __name__ == "__main__":
    main()
