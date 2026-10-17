"""sigmagraph: valuation-weighted sigma-graphs of integers.

Companion code to the notes "Odd Perfect Numbers as Edge-Weighted Digraphs"
(see notes/ in the repository).
"""

from .arith import (abundancy, big_omega, classify, cyclotomic_layers,
                    edge_type, factor, local_sigma, omega, order,
                    primitive_prime_divisors, sigma, sigma0, sigma_r,
                    spoof_sigma, valuation)
from .graph import SigmaGraph
from .divisors import (deficient_twin, divisor_graph_edges,
                       edge_weighted_abundancy, exponent_pattern,
                       hasse_graph_edges, partition_function,
                       top_weighted_degree)
from .data import MULTIPERFECT, DESCARTES_SPOOF

__all__ = [
    "SigmaGraph", "abundancy", "big_omega", "classify", "cyclotomic_layers",
    "edge_type", "factor", "local_sigma", "omega", "order",
    "primitive_prime_divisors", "sigma", "sigma0", "sigma_r", "spoof_sigma",
    "valuation", "deficient_twin", "divisor_graph_edges",
    "edge_weighted_abundancy", "exponent_pattern", "hasse_graph_edges",
    "partition_function", "top_weighted_degree", "MULTIPERFECT",
    "DESCARTES_SPOOF",
]

__version__ = "0.1.0"
