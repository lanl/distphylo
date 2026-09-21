#!/usr/bin/env python3
"""Prune-Deep Forest reconstruction from a precomputed distance matrix.

User command examples:
    python prune_deep_forest_distance.py distances.tsv --outdir prune_deep_output
    python prune_deep_forest_distance.py distances.tsv --outdir prune_deep_output --tau-refine-rounds 1
    python prune_deep_forest_distance.py distances.tsv \
      --outdir prune_deep_output \
      --M "[1.5,2]" \
      --m "[0.2,0.4,0.6]" \
      --tau "[0.005,0.02]"

Input
-----
One full, labeled, symmetric pairwise distance matrix in CSV, TSV, or
whitespace-delimited text. The first row contains column taxon labels and each
subsequent row begins with its row taxon label. For example:

    taxon  A     B     C
    A      0     0.10  0.20
    B      0.10  0     0.15
    C      0.20  0.15  0

The program does not estimate distances from sequences and does not assume a
substitution model. Any model choice used to create the input distance matrix
is external to this program.

Outputs
-------
The most-resolved split-compatible result is always written as the
primary Newick file. If that primary result is a single tree and at least one
split-compatible forest exists, the most-resolved forest is also written as a
secondary ``forest_less_preferred_*.nwk`` file. The grid-search log records
every tested (M,m,tau) point and compares the primary result with the secondary
forest when one is emitted.

Automatic parameter search based on the NJ tree inferred from the supplied
pairwise distance matrix. User can also provide own parameters to use.
--------------------------
* tau center = median positive NJ internal branch length / 4.
  Initial tau multipliers: [0.25, 0.5, 1, 1.5, 2].
* m center = NJ-estimated chord depth using the supplied distance matrix.
  m multipliers: [0.5, 1, 1.5].
* For each (m,tau), M_min = 2*m + 3*tau.
  M multipliers: [1.050, 1.250, 1.500].

Selection
---------
Choose the split-compatible result with the largest total number of nontrivial
splits. For a forest, this is the sum of split counts over all components. A
single tree and a forest compete on the same resolution measure. Ties prefer
fewer components, then earlier grid order.

If the global winner is a single tree, also identify the most-resolved
split-compatible forest (if any) and save it as a secondary answer.

If the selected automatic-grid result uses the smallest tau tested, the search
adds an even smaller tau (half the current minimum) and repeats, up to the
configured refinement limit.
"""

from __future__ import annotations

import argparse
import csv
import io
import math
import re
import time
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Iterable, Sequence

import networkx as nx
import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

EPS = 1e-12


def _split_matrix_line(line: str, delimiter: str | None) -> list[str]:
    """Split one distance-matrix line while preserving CSV/TSV fields."""
    if delimiter is None:
        return line.split()
    return next(csv.reader([line], delimiter=delimiter, skipinitialspace=True))


def read_distance_matrix(path: Path) -> tuple[list[str], np.ndarray, DistanceMatrix]:
    """Read and validate a full labeled pairwise distance matrix.

    Accepted text formats are comma-separated, tab-separated, or whitespace-
    separated. The first row must contain column taxon labels, optionally
    preceded by a corner label such as ``taxon``. Each data row must start with
    its row taxon label and then contain one distance for every column taxon.
    """
    raw_lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if len(raw_lines) < 4:
        raise ValueError("Distance matrix must contain a header and at least 3 taxa")

    first = raw_lines[0]
    delimiter = "," if "," in first else ("\t" if "\t" in first else None)
    rows = [_split_matrix_line(line, delimiter) for line in raw_lines]

    header = [field.strip() for field in rows[0]]
    body = rows[1:]

    # Two supported header styles:
    #   taxon A B C     (explicit corner field)
    #         A B C     (blank CSV/TSV corner field)
    # or, for whitespace-delimited files, simply:
    #   A B C
    if len(header) == len(body) + 1:
        labels = header[1:]
    elif len(header) == len(body):
        labels = header
    else:
        raise ValueError(
            "Distance-matrix header must contain one column label per taxon "
            "(optionally preceded by a corner label such as 'taxon')"
        )

    labels = [label.strip() for label in labels]
    n = len(labels)
    if n < 3:
        raise ValueError("At least 3 taxa are required")
    if any(not label for label in labels):
        raise ValueError("Taxon labels must be non-empty")
    if len(set(labels)) != n:
        raise ValueError("Column taxon labels must be unique")

    row_values: dict[str, list[float]] = {}
    for row_num, row in enumerate(body, start=2):
        if len(row) != n + 1:
            raise ValueError(
                f"Row {row_num} has {len(row) - 1} distances; expected {n}"
            )
        row_label = row[0].strip()
        if row_label in row_values:
            raise ValueError(f"Duplicate row taxon label: {row_label!r}")
        try:
            values = [float(x) for x in row[1:]]
        except ValueError as exc:
            raise ValueError(f"Non-numeric distance in row {row_num} ({row_label!r})") from exc
        row_values[row_label] = values

    if set(row_values) != set(labels):
        missing = sorted(set(labels) - set(row_values))
        extra = sorted(set(row_values) - set(labels))
        raise ValueError(
            f"Row/column taxon labels differ; missing rows={missing}, extra rows={extra}"
        )

    full = np.asarray([row_values[label] for label in labels], dtype=float)
    if full.shape != (n, n):
        raise ValueError(f"Distance matrix must be square; got shape {full.shape}")
    if not np.all(np.isfinite(full)):
        raise ValueError("Distance matrix contains NaN or infinite values")
    if np.any(full < -EPS):
        raise ValueError("Distance matrix contains negative distances")
    if not np.allclose(np.diag(full), 0.0, rtol=0.0, atol=1e-10):
        raise ValueError("Distance-matrix diagonal must be zero")
    if not np.allclose(full, full.T, rtol=1e-8, atol=1e-10):
        max_asym = float(np.max(np.abs(full - full.T)))
        raise ValueError(
            f"Distance matrix must be symmetric; maximum asymmetry is {max_asym:.6g}"
        )

    # Remove tiny floating-point asymmetry and tiny negative roundoff.
    full = 0.5 * (full + full.T)
    full[np.abs(full) < EPS] = 0.0
    np.fill_diagonal(full, 0.0)

    lower = [[float(full[i, j]) for j in range(i + 1)] for i in range(n)]
    dm = DistanceMatrix(names=labels, matrix=lower)
    return labels, full, dm

def infer_nj(dm: DistanceMatrix):
    return DistanceTreeConstructor().nj(dm)


def tree_to_newick(tree) -> str:
    buf = io.StringIO()
    Phylo.write(tree, buf, "newick")
    return buf.getvalue().strip()


def positive_nj_internal_lengths(tree) -> list[float]:
    values = []
    root = tree.root
    for clade in tree.get_nonterminals(order="preorder"):
        if clade is root:
            continue
        length = clade.branch_length
        if length is not None and math.isfinite(length) and length > EPS:
            values.append(float(length))
    return sorted(values)


def pairwise_positive_distances(matrix: np.ndarray) -> np.ndarray:
    vals = matrix[np.triu_indices_from(matrix, k=1)]
    vals = vals[np.isfinite(vals) & (vals > EPS)]
    if vals.size == 0:
        raise ValueError("All pairwise distances are zero; no phylogenetic signal is available")
    return vals


def fmt_num(x: float) -> str:
    """Compact precision for filenames and console output."""
    return f"{x:.6g}"


def fmt_log(x: float) -> str:
    """Readable fixed precision for the grid-search log only."""
    return f"{x:.3f}"


def newick_for_log(newick: str) -> str:
    """Round branch lengths in a Newick string to 3 decimals for log display."""
    def repl(match):
        return ":" + f"{float(match.group(1)):.3f}"
    return re.sub(r":([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)", repl, newick)


def theorem_valid(M: float, m: float, tau: float) -> bool:
    return M > 3.0 * tau and m > 3.0 * tau and M > 2.0 * m + 3.0 * tau


def nj_chord_depth(
    nj_tree,
    labels: Sequence[str],
    distance_matrix: np.ndarray,
) -> tuple[float, list[float]]:
    """Estimate chord depth from NJ internal splits and supplied distances.

    For each unique nontrivial split induced by an NJ internal edge, compute the
    minimum supplied leaf-to-leaf distance crossing that split. The NJ
    chord-depth proxy is the maximum of those split-specific minimum distances.
    """
    label_to_idx = {label: i for i, label in enumerate(labels)}
    all_labels = set(labels)
    seen: set[tuple[str, ...]] = set()
    depths: list[float] = []

    for clade in nj_tree.find_clades(order="preorder"):
        if clade is nj_tree.root:
            continue
        side = {t.name for t in clade.get_terminals() if t.name is not None}
        other = all_labels - side
        if len(side) < 2 or len(other) < 2:
            continue

        side_key = tuple(sorted(side))
        other_key = tuple(sorted(other))
        key = side_key if (len(side_key), side_key) <= (len(other_key), other_key) else other_key
        if key in seen:
            continue
        seen.add(key)

        crossing = [
            distance_matrix[label_to_idx[a], label_to_idx[b]]
            for a in side
            for b in other
        ]
        if crossing:
            depths.append(float(min(crossing)))

    if depths:
        depths.sort()
        return float(max(depths)), depths

    # Very small/star-like NJ trees may not have a nontrivial internal split.
    # Use the largest positive pairwise distance only as a fallback scale.
    pairwise = pairwise_positive_distances(distance_matrix)
    fallback = float(np.max(pairwise))
    return fallback, [fallback]


def build_auto_grid(
    tau_values: Sequence[float],
    chord_depth: float,
) -> list[tuple[float, float, float]]:
    """Create the automatic grid and discard non-Theorem-1 combinations."""
    m_values = [chord_depth * f for f in (0.5, 1.0, 1.5)]
    triples: list[tuple[float, float, float]] = []
    for tau in tau_values:
        for m in m_values:
            M_min = 2.0 * m + 3.0 * tau
            for mult in (1.050, 1.250, 1.500):
                M = float(mult * M_min)
                if theorem_valid(M, m, tau):
                    triples.append((M, float(m), float(tau)))
    return triples


def auto_grid(
    nj_tree,
    labels: Sequence[str],
    distance_matrix: np.ndarray,
) -> tuple[list[tuple[float, float, float]], dict]:
    """Build the initial NJ-informed grid."""
    pairwise = pairwise_positive_distances(distance_matrix)
    internal = positive_nj_internal_lengths(nj_tree)

    if internal:
        tau_center = float(np.median(internal)) / 4.0
        tau_source = "median positive NJ internal branch / 4"
    else:
        tau_center = float(np.quantile(pairwise, 0.10)) / 4.0
        tau_source = "10th percentile positive pairwise distance / 4 (NJ had no positive internal edge)"
    tau_center = max(tau_center, 1e-12)

    tau_multipliers = (0.25, 0.5, 1.0, 1.5, 2.0)
    tau_values = sorted({tau_center * f for f in tau_multipliers})

    chord_depth, chord_depths = nj_chord_depth(nj_tree, labels, distance_matrix)
    m_multipliers = (0.5, 1.0, 1.5)
    m_values = [chord_depth * f for f in m_multipliers]
    M_multipliers = (1.050, 1.250, 1.500)

    raw_grid_points = len(tau_values) * len(m_values) * len(M_multipliers)
    grid = build_auto_grid(tau_values, chord_depth)
    metadata = {
        "tau_source": tau_source,
        "tau_center": tau_center,
        "tau_multipliers": list(tau_multipliers),
        "tau_values_initial": list(tau_values),
        "tau_values_final": list(tau_values),
        "tau_refinement_rounds_used": 0,
        "nj_internal_lengths": internal,
        "nj_chord_depth": chord_depth,
        "nj_split_chord_depths": chord_depths,
        "m_multipliers": list(m_multipliers),
        "m_values": m_values,
        "M_multipliers": list(M_multipliers),
        "grid_points_initial_raw": raw_grid_points,
        "grid_points_initial": len(grid),
        "grid_points_initial_discarded": raw_grid_points - len(grid),
        "grid_points_raw_total": raw_grid_points,
        "grid_points_discarded_total": raw_grid_points - len(grid),
    }
    return grid, metadata


def manual_grid(
    Ms: Sequence[float],
    ms: Sequence[float],
    taus: Sequence[float],
) -> list[tuple[float, float, float]]:
    """Custom (manual input) grid restricted to strict Theorem-1 combinations."""
    return [
        (float(M), float(m), float(tau))
        for M, m, tau in product(Ms, ms, taus)
        if theorem_valid(float(M), float(m), float(tau))
    ]

# FOREST construction functions
def construct_clustering_graph(distance_matrix: np.ndarray, m: float) -> nx.Graph:
    n = distance_matrix.shape[0]
    G = nx.Graph()
    G.add_nodes_from(range(n))
    for i in range(n):
        for j in range(i + 1, n):
            if distance_matrix[i, j] < m:
                G.add_edge(i, j)
    return G


def calculate_phi(u: int, v: int, w: int, distance_matrix: np.ndarray) -> float:
    return 0.5 * (
        distance_matrix[u, v] + distance_matrix[u, w] - distance_matrix[v, w]
    )


def mini_contractor(
    component_nodes: set[int],
    distance_matrix: np.ndarray,
    leaves: tuple[int, int],
    M: float,
    tau: float,
) -> list[tuple[set[int], set[int]]]:
    u, v = leaves
    B = {
        w
        for w in component_nodes
        if max(distance_matrix[u, w], distance_matrix[v, w]) < M
    }
    B.update((u, v))
    S = B - {u}
    bipartitions: list[tuple[set[int], set[int]]] = []
    x_prev = u
    phi = {w: calculate_phi(u, v, w, distance_matrix) for w in S}
    phi[u] = 0.0

    while S:
        x0 = min(S, key=lambda w: phi[w])
        if phi[x0] - phi[x_prev] >= 2.0 * tau:
            bipartitions.append((B - S.copy(), S.copy()))
        S.remove(x0)
        x_prev = x0
    return bipartitions


def extender(
    graph: nx.Graph,
    bipartitions: Iterable[tuple[set[int], set[int]]],
) -> list[tuple[set[int], set[int]]]:
    extended = []
    for psi_u0, psi_v0 in bipartitions:
        psi_u, psi_v = set(psi_u0), set(psi_v0)
        K = graph.copy()
        K.remove_edges_from(
            (a, b) for a in psi_u for b in psi_v if K.has_edge(a, b)
        )

        if any(nx.has_path(K, a, b) for a in psi_u for b in psi_v):
            continue

        out_u, out_v = set(psi_u), set(psi_v)
        valid = True
        for w in set(K.nodes) - (psi_u | psi_v):
            to_u = any(nx.has_path(K, w, a) for a in psi_u)
            to_v = any(nx.has_path(K, w, b) for b in psi_v)
            if to_u and not to_v:
                out_u.add(w)
            elif to_v and not to_u:
                out_v.add(w)
            else:
                valid = False
                break
        if valid:
            extended.append((out_u, out_v))
    return extended

# After bipartitions constructed, look at splits and compatibility
def canonical_split(split: tuple[set[int], set[int]]) -> tuple[frozenset[int], frozenset[int]]:
    A, B = map(frozenset, split)
    ka = (len(A), tuple(sorted(A)))
    kb = (len(B), tuple(sorted(B)))
    return (A, B) if ka <= kb else (B, A)


def split_is_nontrivial(split: tuple[Iterable[int], Iterable[int]]) -> bool:
    A, B = split
    return len(A) >= 2 and len(B) >= 2


def splits_compatible(a, b) -> bool:
    A, B = map(set, a)
    C, D = map(set, b)
    return not (A & C and A & D and B & C and B & D)


def all_compatible(splits: Sequence[tuple[frozenset[int], frozenset[int]]]) -> bool:
    for i in range(len(splits)):
        for j in range(i + 1, len(splits)):
            if not splits_compatible(splits[i], splits[j]):
                return False
    return True


def quote_newick_label(label: str) -> str:
    if re.fullmatch(r"[A-Za-z0-9_.-]+", label):
        return label
    return "'" + label.replace("'", "''") + "'"


def component_newick(
    component: set[int],
    splits: Sequence[tuple[frozenset[int], frozenset[int]]],
    labels: Sequence[str],
) -> str:
    taxa = set(component)
    if len(taxa) == 1:
        return quote_newick_label(labels[next(iter(taxa))]) + ";"
    if len(taxa) == 2:
        return "(" + ",".join(quote_newick_label(labels[i]) for i in sorted(taxa)) + ");"

    ref = min(taxa)
    clusters: set[frozenset[int]] = set()
    for A0, B0 in splits:
        A, B = set(A0), set(B0)
        if A | B != taxa or A & B:
            continue
        side = B if ref in A else A
        if 1 < len(side) < len(taxa):
            clusters.add(frozenset(side))

    cluster_list = sorted(clusters, key=lambda c: (len(c), tuple(sorted(c))))

    def render_cluster(C: frozenset[int]) -> str:
        proper = [D for D in cluster_list if D < C]
        maximal = [D for D in proper if not any(D < E < C for E in proper)]
        covered = set().union(*(set(D) for D in maximal)) if maximal else set()
        parts = [render_cluster(D) for D in sorted(maximal, key=lambda x: min(x))]
        parts.extend(quote_newick_label(labels[i]) for i in sorted(set(C) - covered))
        return "(" + ",".join(parts) + ")"

    top = [C for C in cluster_list if not any(C < D for D in cluster_list)]
    covered_top = set().union(*(set(C) for C in top)) if top else set()
    parts = [render_cluster(C) for C in sorted(top, key=lambda x: min(x))]
    parts.extend(quote_newick_label(labels[i]) for i in sorted(taxa - covered_top))
    return "(" + ",".join(parts) + ");"


@dataclass
class Candidate:
    index: int
    M: float
    m: float
    tau: float
    split_compatible: bool
    num_components: int
    component_sizes: list[int]
    split_count: int
    component_split_counts: list[int]
    forest_newicks: list[str]

    @property
    def selection_key(self):
        # One tree and a forest compete on the same resolution metric.
        # For a forest, split_count is already the sum over components.
        # All evaluated points already satisfy Theorem 1.
        # On an exact resolution tie, prefer fewer components, then earlier grid order.
        return (
            self.split_count,
            -self.num_components,
            -self.index,
        )


def evaluate_grid_point(
    index: int,
    M: float,
    m: float,
    tau: float,
    labels: Sequence[str],
    distances: np.ndarray,
) -> Candidate:
    G = construct_clustering_graph(distances, m)
    components = [set(c) for c in nx.connected_components(G)]
    components.sort(key=lambda c: min(c))
    forest_newicks: list[str] = []
    component_split_counts: list[int] = []

    for component in components:
        if len(component) <= 2:
            forest_newicks.append(component_newick(component, [], labels))
            component_split_counts.append(0)
            continue

        subgraph = G.subgraph(component).copy()
        unique: set[tuple[frozenset[int], frozenset[int]]] = set()

        # Prune-Deep Mini Contractor runs only on m-clustering graph edges.
        for leaves in subgraph.edges():
            mini = mini_contractor(component, distances, leaves, M, tau)
            for split in mini:
                A, B = split
                if A | B == component and not (A & B):
                    candidates = [split]
                else:
                    candidates = extender(subgraph, [split])
                for candidate in candidates:
                    A2, B2 = candidate
                    if A2 | B2 == component and not (A2 & B2):
                        unique.add(canonical_split(candidate))

        nontrivial = sorted(
            [s for s in unique if split_is_nontrivial(s)],
            key=lambda s: (len(s[0]), tuple(sorted(s[0])), tuple(sorted(s[1]))),
        )
        if not all_compatible(nontrivial):
            return Candidate(
                index=index, M=M, m=m, tau=tau,
                split_compatible=False,
                num_components=len(components),
                component_sizes=[len(c) for c in components],
                split_count=0, component_split_counts=[], forest_newicks=[]
            )

        forest_newicks.append(component_newick(component, nontrivial, labels))
        component_split_counts.append(len(nontrivial))

    return Candidate(
        index=index,
        M=M,
        m=m,
        tau=tau,
        split_compatible=True,
        num_components=len(components),
        component_sizes=[len(c) for c in components],
        split_count=sum(component_split_counts),
        component_split_counts=component_split_counts,
        forest_newicks=forest_newicks,
    )


def parameter_filename(M: float, m: float, tau: float) -> str:
    return f"main_M{fmt_num(M)}m{fmt_num(m)}tau{fmt_num(tau)}.nwk"


def secondary_forest_filename(M: float, m: float, tau: float) -> str:
    return f"forest_less_preferred_M{fmt_num(M)}m{fmt_num(m)}tau{fmt_num(tau)}.nwk"


def write_log(
    path: Path,
    input_path: Path,
    labels: Sequence[str],
    nj_newick: str,
    grid_meta: dict,
    candidates: Sequence[Candidate],
    best: Candidate,
    output_path: Path,
    secondary_forest: Candidate | None,
    secondary_forest_path: Path | None,
    elapsed: float,
    manual: bool,
) -> None:
    with path.open("w", encoding="utf-8") as fh:
        fh.write("Prune-Deep Forest distance-matrix grid search\n")
        fh.write(f"input_distance_matrix\t{input_path}\n")
        fh.write(f"taxa\t{len(labels)}\n")
        fh.write("distance_source\tuser_supplied_precomputed_matrix\n")
        fh.write("distance_model\tunspecified_by_program\n")
        fh.write(f"grid_mode\t{'custom' if manual else 'NJ-informed automatic'}\n")

        if manual:
            fh.write("custom_M_values\t" + ",".join(fmt_log(x) for x in grid_meta["Ms"]) + "\n")
            fh.write("custom_m_values\t" + ",".join(fmt_log(x) for x in grid_meta["ms"]) + "\n")
            fh.write("custom_tau_values\t" + ",".join(fmt_log(x) for x in grid_meta["taus"]) + "\n")
        else:
            fh.write(
                "NJ_inferred_internal_branches\t"
                + ",".join(fmt_log(x) for x in grid_meta["nj_internal_lengths"])
                + "\n"
            )
            fh.write(
                f"tau\tmedian(NJ internal branches)/4={fmt_log(grid_meta['tau_center'])}; "
                "multipliers=[0.250,0.500,1.000,1.500,2.000]\n"
            )
            fh.write(
                "tau_values_initial\t"
                + ",".join(fmt_log(x) for x in grid_meta["tau_values_initial"])
                + "\n"
            )
            fh.write(
                "tau_values_final\t"
                + ",".join(fmt_log(x) for x in grid_meta["tau_values_final"])
                + "\n"
            )
            fh.write(
                f"tau_refinement_rounds_used\t{grid_meta['tau_refinement_rounds_used']}\n"
            )
            fh.write(
                f"m\tNJ chord depth={fmt_log(grid_meta['nj_chord_depth'])}; "
                "multipliers=[0.500,1.000,1.500]\n"
            )
            fh.write(
                "NJ_split_chord_depths\t"
                + ",".join(fmt_log(x) for x in grid_meta["nj_split_chord_depths"])
                + "\n"
            )
            fh.write(
                "m_values\t"
                + ",".join(fmt_log(x) for x in grid_meta["m_values"])
                + "\n"
            )
            fh.write(
                "M\tM_base=2*m+3*tau; multipliers=[1.050, 1.250, 1.500]\n"
            )
            fh.write(f"grid_points_initial_raw\t{grid_meta['grid_points_initial_raw']}\n")
            fh.write(f"grid_points_initial_discarded\t{grid_meta['grid_points_initial_discarded']}\n")
            fh.write(f"grid_points_initial_tested\t{grid_meta['grid_points_initial']}\n")

        fh.write(f"nj_tree\t{newick_for_log(nj_newick)}\n")
        fh.write("theorem_condition\tM>3*tau; m>3*tau; M>2*m+3*tau\n")
        fh.write(
            "selection\tmaximum total compatible nontrivial split_count; "
            "forest split_count is summed across components; tie -> fewer components, then earlier grid point\n"
        )
        if not manual:
            fh.write(f"grid_points_raw_total\t{grid_meta['grid_points_raw_total']}\n")
            fh.write(f"grid_points_discarded_by_theorem_filter\t{grid_meta['grid_points_discarded_total']}\n")
        else:
            fh.write(f"grid_points_raw_total\t{grid_meta['raw_count']}\n")
            fh.write(f"grid_points_discarded_by_theorem_filter\t{grid_meta['discarded_count']}\n")
        fh.write(f"grid_points_tested\t{len(candidates)}\n")
        fh.write(f"elapsed_seconds\t{elapsed:.3f}\n\n")

        fh.write(
            "grid_index\tM\tm\ttau\tsplit_compatible\t"
            "num_components\tcomponent_sizes\ttotal_split_count\tcomponent_split_counts\n"
        )
        for c in candidates:
            fh.write(
                f"{c.index}\t{fmt_log(c.M)}\t{fmt_log(c.m)}\t{fmt_log(c.tau)}\t"
                f"{int(c.split_compatible)}\t{c.num_components}\t"
                f"{','.join(map(str, c.component_sizes))}\t{c.split_count}\t"
                f"{','.join(map(str, c.component_split_counts))}\n"
            )

        fh.write("\nSELECTED\n")
        fh.write(f"M\t{fmt_log(best.M)}\n")
        fh.write(f"m\t{fmt_log(best.m)}\n")
        fh.write(f"tau\t{fmt_log(best.tau)}\n")
        fh.write(f"num_components\t{best.num_components}\n")
        fh.write("component_sizes\t" + ",".join(map(str, best.component_sizes)) + "\n")
        fh.write(f"total_split_count\t{best.split_count}\n")
        fh.write(
            "component_split_counts\t"
            + ",".join(map(str, best.component_split_counts))
            + "\n"
        )
        fh.write(f"newick_file\t{output_path}\n")
        fh.write(f"selected_result_type\t{'forest' if best.num_components > 1 else 'single_tree'}\n")

        if secondary_forest is not None and secondary_forest_path is not None:
            split_difference = best.split_count - secondary_forest.split_count
            fh.write("\nSECONDARY_MOST_RESOLVED_FOREST\n")
            fh.write("note\tPrimary answer is a single tree; this is the most-resolved compatible forest found.\n")
            if split_difference > 0:
                fh.write(
                    f"resolution_note\tForest is less resolved than selected single tree by {split_difference} nontrivial split(s).\n"
                )
            else:
                fh.write(
                    "resolution_note\tForest has the same total split count as the selected single tree but lost the tie-break.\n"
                )
            fh.write(f"M\t{fmt_log(secondary_forest.M)}\n")
            fh.write(f"m\t{fmt_log(secondary_forest.m)}\n")
            fh.write(f"tau\t{fmt_log(secondary_forest.tau)}\n")
            fh.write(f"num_components\t{secondary_forest.num_components}\n")
            fh.write(
                "component_sizes\t"
                + ",".join(map(str, secondary_forest.component_sizes))
                + "\n"
            )
            fh.write(f"total_split_count\t{secondary_forest.split_count}\n")
            fh.write(
                "component_split_counts\t"
                + ",".join(map(str, secondary_forest.component_split_counts))
                + "\n"
            )
            fh.write(f"newick_file\t{secondary_forest_path}\n")


def parse_float_list(tokens: Sequence[str], option_name: str) -> list[float]:
    """Accept manual grid parameters as --m 0.1 0.2, --m 0.1,0.2, or --m '[0.1,0.2]'."""
    text = " ".join(tokens).strip().replace("[", " ").replace("]", " ")
    pieces = [p for p in re.split(r"[\s,]+", text) if p]
    try:
        values = [float(p) for p in pieces]
    except ValueError as exc:
        raise SystemExit(f"ERROR: {option_name} must contain only numeric values") from exc
    if not values:
        raise SystemExit(f"ERROR: {option_name} requires at least one numeric value")
    return values


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Reconstruct one most-resolved Prune-Deep tree/forest from a precomputed "
            "pairwise distance matrix. Without custom grid options, M/m/tau are searched "
            "automatically from NJ scales."
        )
    )
    p.add_argument(
        "distance_matrix", type=Path,
        help="Full labeled symmetric distance matrix (CSV, TSV, or whitespace-delimited text)",
    )
    p.add_argument("--outdir", type=Path, default=Path("./prune_deep_output"), help="Output directory")
    p.add_argument("--output", type=str, help="Newick output basename. Default: M{M}m{m}tau{tau}.nwk")
    p.add_argument("--log", type=str, default="grid_search.log", help="Grid-search log basename")
    p.add_argument(
        "--M", dest="Ms_raw", type=str, nargs="+",
        help="Custom M list, e.g. --M 0.5 1 2 or --M '[0.5,1,2]'; requires --m and --tau",
    )
    p.add_argument(
        "--m", dest="ms_raw", type=str, nargs="+",
        help="Custom m list, e.g. --m 0.1 0.2 or --m '[0.1,0.2]'; requires --M and --tau",
    )
    p.add_argument(
        "--tau", dest="taus_raw", type=str, nargs="+",
        help="Custom tau list, e.g. --tau 0.01 0.02 or --tau '[0.01,0.02]'; requires --M and --m",
    )
    p.add_argument(
        "--tau-refine-rounds", type=int, default=3,
        help=(
            "Automatic mode only: if the selected result uses the smallest tau, "
            "add tau/2 and rerun, up to this many rounds (default: 3)"
        ),
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    distance_matrix_path = args.distance_matrix.expanduser().resolve()
    if not distance_matrix_path.is_file():
        raise SystemExit(f"ERROR: distance matrix file not found: {distance_matrix_path}")
    if args.tau_refine_rounds < 0:
        raise SystemExit("ERROR: --tau-refine-rounds must be >= 0")

    custom_flags = [args.Ms_raw is not None, args.ms_raw is not None, args.taus_raw is not None]
    if any(custom_flags) and not all(custom_flags):
        raise SystemExit("ERROR: custom grid requires all three: --M, --m, and --tau")

    start = time.perf_counter()
    try:
        labels, distances, dm = read_distance_matrix(distance_matrix_path)
    except ValueError as exc:
        raise SystemExit(f"ERROR: {exc}") from exc
    nj_tree = infer_nj(dm)
    nj_newick = tree_to_newick(nj_tree)

    manual = all(custom_flags)
    if manual:
        Ms = parse_float_list(args.Ms_raw, "--M")
        ms = parse_float_list(args.ms_raw, "--m")
        taus = parse_float_list(args.taus_raw, "--tau")
        raw_count = len(Ms) * len(ms) * len(taus)
        grid = manual_grid(Ms, ms, taus)
        if not grid:
            raise SystemExit(
                "ERROR: no custom grid point satisfies Theorem 1: "
                "M>3*tau, m>3*tau, and M>2*m+3*tau"
            )
        grid_meta = {
            "Ms": Ms, "ms": ms, "taus": taus,
            "raw_count": raw_count,
            "discarded_count": raw_count - len(grid),
        }
    else:
        grid, grid_meta = auto_grid(nj_tree, labels, distances)

    candidates: list[Candidate] = []
    seen_grid: set[tuple[float, float, float]] = set()

    def evaluate_new_grid(points: Sequence[tuple[float, float, float]]) -> None:
        for M, m, tau in points:
            key = (round(M, 14), round(m, 14), round(tau, 14))
            if key in seen_grid:
                continue
            seen_grid.add(key)
            candidates.append(
                evaluate_grid_point(
                    len(candidates) + 1, M, m, tau, labels, distances
                )
            )

    evaluate_new_grid(grid)

    def compatible_candidates() -> list[Candidate]:
        return [c for c in candidates if c.split_compatible]

    compatible = compatible_candidates()
    if not compatible:
        raise SystemExit("ERROR: no grid point produced a mutually compatible split set")
    best = max(compatible, key=lambda c: c.selection_key)

    # Adaptive lower-tau refinement: if the winner is at the lower edge of the
    # automatic tau grid, search another factor of two lower. Repeat up to the
    # user-configured limit.
    if not manual:
        tau_values = list(grid_meta["tau_values_final"])
        rounds_used = 0
        while rounds_used < args.tau_refine_rounds:
            smallest_tau = min(tau_values)
            if not math.isclose(best.tau, smallest_tau, rel_tol=1e-10, abs_tol=1e-15):
                break
            new_tau = smallest_tau / 2.0
            if new_tau <= 0.0:
                break
            tau_values.append(new_tau)
            tau_values = sorted(set(tau_values))
            refine_raw_count = len(grid_meta["m_multipliers"]) * len(grid_meta["M_multipliers"])
            refine_grid = build_auto_grid([new_tau], grid_meta["nj_chord_depth"])
            grid_meta["grid_points_raw_total"] += refine_raw_count
            grid_meta["grid_points_discarded_total"] += refine_raw_count - len(refine_grid)
            evaluate_new_grid(refine_grid)
            rounds_used += 1
            compatible = compatible_candidates()
            best = max(compatible, key=lambda c: c.selection_key)

        grid_meta["tau_values_final"] = tau_values
        grid_meta["tau_refinement_rounds_used"] = rounds_used

    outdir = args.outdir.expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    output_name = args.output or parameter_filename(best.M, best.m, best.tau)
    output_path = outdir / output_name
    log_path = outdir / args.log

    with output_path.open("w", encoding="utf-8") as fh:
        for newick in best.forest_newicks:
            fh.write(newick.rstrip() + "\n")

    secondary_forest: Candidate | None = None
    secondary_forest_path: Path | None = None
    if best.num_components == 1:
        forest_candidates = [
            c for c in compatible if c.num_components > 1
        ]
        if forest_candidates:
            secondary_forest = max(forest_candidates, key=lambda c: c.selection_key)
            secondary_forest_path = outdir / secondary_forest_filename(
                secondary_forest.M, secondary_forest.m, secondary_forest.tau
            )
            with secondary_forest_path.open("w", encoding="utf-8") as fh:
                for newick in secondary_forest.forest_newicks:
                    fh.write(newick.rstrip() + "\n")

    elapsed = time.perf_counter() - start
    write_log(
        log_path, distance_matrix_path, labels, nj_newick,
        grid_meta, candidates, best, output_path, secondary_forest,
        secondary_forest_path, elapsed, manual,
    )

    print(f"Selected M={fmt_num(best.M)} m={fmt_num(best.m)} tau={fmt_num(best.tau)}")
    print(
        f"Components={best.num_components} sizes={best.component_sizes} "
        f"total_nontrivial_splits={best.split_count}"
    )
    if not manual:
        print(
            f"Tau refinement rounds={grid_meta['tau_refinement_rounds_used']} "
            f"final_min_tau={fmt_num(min(grid_meta['tau_values_final']))}"
        )
    print(f"Newick: {output_path}")
    if secondary_forest is not None and secondary_forest_path is not None:
        print(
            f"Secondary forest: {secondary_forest_path} "
            f"(total_nontrivial_splits={secondary_forest.split_count}, "
            f"primary={best.split_count})"
        )
    print(f"Log:    {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
