"""LaTeX and TikZ output for sigma-graphs."""

from __future__ import annotations

import math

from .graph import SigmaGraph


def latex_table(G: SigmaGraph) -> str:
    """Vertices with exponent a, in-degree b, out-degree, and defect d."""
    d = G.defect()
    rows = [r"\begin{tabular}{@{}rrrrr@{}}", r"\toprule",
            r"prime & $a$ & in-degree $b$ & out-degree & defect $d$\\",
            r"\midrule"]
    for v in G.vertices:
        rows.append(f"${v}$ & ${G.a[v]}$ & ${G.in_degree(v)}$ & "
                    f"${G.out_degree(v)}$ & ${d[v]}$\\\\")
    rows += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(rows)


def tikz(G: SigmaGraph, radius: float = 2.2) -> str:
    """TikZ picture of G on a circle; boundary vertices are dashed."""
    V = G.vertices
    lines = [r"\begin{tikzpicture}[>=Stealth,",
             r"  v/.style={circle,draw,minimum size=8mm,inner sep=0pt},",
             r"  bv/.style={circle,draw,dashed,minimum size=8mm,inner sep=0pt}]"]
    m = len(V)
    for i, v in enumerate(V):
        ang = 90 + 360 * i / m
        x, y = radius * math.cos(math.radians(ang)), radius * math.sin(math.radians(ang))
        style = "bv" if G.a[v] == 0 else "v"
        lines.append(rf"\node[{style}] (p{v}) at ({x:.2f},{y:.2f}) {{${v}$}};")
    for s, r, e in G.edges():
        bend = ",bend left=15" if G.w(r, s) else ""
        lines.append(rf"\draw[->{bend}] (p{s}) to node[auto,font=\small]{{{e}}} (p{r});")
    lines.append(r"\end{tikzpicture}")
    return "\n".join(lines)
