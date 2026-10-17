# sigmagraph — σ-graphs of integers and odd perfect numbers

Companion code and notes for **"Odd Perfect Numbers as Edge-Weighted
Digraphs"** (David Joyner, working notes prepared with Claude, Gemini and
ChatGPT). The notes are in [`notes/`](notes/), as LaTeX source and PDF.

The **σ-graph** of an integer `n > 1` is a weighted digraph on primes. Its
vertices are the primes dividing `n·σ(n)`. The edge `s → r` has weight
`v_r(σ(s^a))`, where `s^a` is the exact power of `s` in `n`. So an edge says
"the prime `r` occurs in the local divisor sum at `s`", and its weight says
how often.

This package computes these graphs exactly and reproduces every numerical
example in the notes as a test.

## What it computes

For any integer `n > 1`:

* the **extended** and **internal** σ-graphs, with boundary vertices
  (primes dividing `σ(n)` but not `n`);
* the **defect vector** `d(n)`, with the check `σ(n)/n = ∏ r^{d_r}`
  (`n` is perfect exactly when `d(n) = e_2`);
* the **Poisson operator** `L = diag(a) − W`, with the checks
  `Lℓ = −η` and `1ᵀL = −d(n)ᵀ`, where `ℓ = (log p)` and `η` holds the local
  log-abundancies;
* the **order/LTE type** of each edge between odd primes (order-one,
  Wieferich-type, multiplier);
* strong connectivity and strongly connected components;
* the canonical **arithmetical structure** `(R, S)` with `S` = in-degrees
  (Lorenzini, extended to digraphs), and its **critical group** via the
  Smith normal form;
* the quadratic-form ("energy") identity;
* divisor graphs `D(n)`, Hasse graphs `H(n)`, their partition function, the
  edge-weighted form of the abundancy, and "deficient twins" with the same
  unlabeled divisor graph;
* LaTeX tables and TikZ pictures.

All arithmetic is exact (integers, rationals, Smith normal form). Floating
point is used only for logarithms in the Poisson and energy checks.

## Install

```bash
git clone https://github.com/<user>/opn-sigma-graphs.git
cd opn-sigma-graphs
python -m pip install -e ".[test]"
```

The only runtime dependency is [SymPy](https://www.sympy.org).

## Use

```bash
$ python -m sigmagraph 120
n = 120 = 2^3 * 3 * 5
sigma(n) = 360,  h(n) = 3  (3-perfect)
vertices: [2, 3, 5]   boundary: []
edges (source -> target : weight):
  2 -> 3 : 1
  2 -> 5 : 1
  3 -> 2 : 2
  5 -> 2 : 1
  5 -> 3 : 1
defect d(n): {2: 0, 3: 1, 5: 0}
strongly connected: True
arithmetical structure R = {2: 1, 3: 1, 5: 2}
                       S = {2: 3, 3: 2, 5: 1}
critical group invariant factors: [2]
```

Other options:

```bash
python -m sigmagraph 945 --json       # everything, as JSON
python -m sigmagraph 945 --latex      # vertex table for LaTeX
python -m sigmagraph 15  --tikz       # TikZ picture (boundary vertices dashed)
python -m sigmagraph 945 --internal   # internal graph only
python examples/survey.py --count 20  # survey of small multiperfect numbers
```

From Python:

```python
from sigmagraph import SigmaGraph, edge_type

G = SigmaGraph(945)
G.edges()          # [(3, 2, 3), (3, 5, 1), (5, 2, 1), (5, 3, 1), (7, 2, 3)]
G.defect()         # {2: 7, 3: -2, 5: 0, 7: -1}
G.check_defect()   # True: 128/63 = 2^7 3^-2 7^-1

edge_type(7, 8, 3) # order-one edge of weight 2: sigma(7^8) = 3^2 * 19 * 37 * 1063
```

## Tests

```bash
pytest -q
```

The tests in [`tests/test_notes.py`](tests/test_notes.py) reproduce the
numerical examples of the notes: divisor functions, the catalog of perfect,
abundant and deficient numbers, the Descartes spoof, the order/LTE
classification (including the counterexamples to Theorem 5.1 of Rodgers'
preprint), the σ-graphs of 6, 12, 15, 28, 120 and 945, the defect table, the
Poisson operators, the arithmetical structures and critical groups of small
multiperfect numbers, and the divisor-graph examples.

## Data

`sigmagraph.MULTIPERFECT` holds the first 26 multiperfect numbers after 1
(the start of OEIS A007691); the tests check each one. For the full list of
known multiperfect numbers, download a published table and read it with
`sigmagraph.data.load_multiperfect(path)`.

## Status and scope

This is exploratory research code. Nothing here proves or claims progress
toward proving that odd perfect numbers do not exist. The notes explain which
statements are classical, which are new but elementary, and which are open.

## License

Code: MIT (see [`LICENSE`](LICENSE)). The notes in [`notes/`](notes/) are
© David Joyner.
