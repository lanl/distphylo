# distphylo
O# (O4745) Inferring trees from limited data

Reconstruct a complete phylogenetic tree or a forest of disjoint trees from a multiple-sequence alignment.

Reference:

The forest algorithm is based on:

> Daskalakis, C., Mossel, E. and Roch, S., 2011. Phylogenies without branch bounds: Contracting the short, pruning the deep. SIAM Journal on Discrete Mathematics, 25(2), pp.872-893.

Input: alignment file, true tree (Newick) for comparison, parameters m, M, and $\tau$. Edit main() function to make changes.

Output: summary tsv files, filename pattern (per replicate):

grid_summary_ntips{ntips}_{replicate}_k{sequence_length}_sorted.tsv

Example command (pass ntips and sequence length, change m, M, $\tau$ in main()):
```
python forest_algorithm.py --n 128 --k 500
```
Output: grid_summary_ntips128_1_k500_sorted.tsv and grid_summary_ntips128_2_k500_sorted.tsv files.

# getforest branch

This branch repository provides two user-facing implementations of the phylogenetic forest reconstruction workflow:

1. **`prune_deep_forest.py`** — starts from an aligned nucleotide FASTA file and estimates pairwise distances with JC69.
2. **`prune_deep_forest_distance.py`** — starts from a given pairwise distance matrix and therefore does not assume a substitution model inside the program.

Both scripts infer a Neighbor Joining (NJ) tree to obtain scales for the automatic parameter search, run Prune-Deep over theorem-valid $(M,m,\tau)$ combinations, and write the most-resolved split-compatible tree or forest. The user can also provide their own parameter values to avoid NJ based values.

---

## Requirements

Python 3.10+ is recommended.

Required Python packages:

```bash
pip install numpy networkx biopython
```
---

### Option 1: aligned FASTA input

Use:

```text
prune_deep_forest.py
```

Pipeline:

```text
aligned FASTA
    -> JC69 pairwise distance matrix
    -> Neighbor Joining tree
    -> automatic or user-specified (M, m, tau) grid
    -> Prune-Deep reconstruction
    -> selected Newick tree/forest + grid-search log
```

The JC69 assumption enters only when converting the alignment into a pairwise evolutionary distance matrix. NJ and the Prune-Deep reconstruction then operate on that distance matrix.

### Option 2: precomputed distance matrix input

Use:

```text
prune_deep_forest_distance.py
```

Pipeline:

```text
precomputed pairwise distance matrix
    -> Neighbor Joining tree
    -> automatic or user-specified (M, m, tau) grid
    -> Prune-Deep reconstruction
    -> selected Newick tree/forest + grid-search log
```

This option does not estimate distances from sequences and does not assume any substitution model internally. Any evolutionary model or distance-estimation method used to create the matrix is external to this program.

---

## 1. Alignment-input version: `prune_deep_forest.py`

### Input

One aligned nucleotide FASTA file with:

- at least 3 sequences;
- equal sequence lengths;
- unique sequence IDs.

This implementation calculates JC69 distances.

### Basic automatic-grid example

```bash
python prune_deep_forest.py alignments/aln_20tips_k100.fa \
    --outdir prune_deep_output_20tips_k100
```

The program calculates JC69 distances, infers an NJ tree, builds the automatic parameter grid, runs Prune-Deep, and selects the most-resolved split-compatible result.

### Limit tau refinement

The default maximum is 3 lower-tau refinement rounds. To allow only one:

```bash
python prune_deep_forest.py alignments/aln_20tips_k100.fa \
    --outdir prune_deep_output_20tips_k100 \
    --tau-refine-rounds 1
```

To disable adaptive lower-tau refinement:

```bash
python prune_deep_forest.py alignments/aln_20tips_k100.fa \
    --outdir prune_deep_output_20tips_k100 \
    --tau-refine-rounds 0
```

### Custom parameter grid

All three options `--M`, `--m`, and `--tau` must be supplied together.

Bracket-list form:

```bash
python prune_deep_forest.py example_alignment.fa \
    --outdir prune_deep_output \
    --M "[1.5,2.0]" \
    --m "[0.2,0.4,0.6]" \
    --tau "[0.005,0.02]"
```

Space-separated form also works:

```bash
python prune_deep_forest.py example_alignment.fa \
    --outdir prune_deep_output \
    --M 1.5 2.0 \
    --m 0.2 0.4 0.6 \
    --tau 0.005 0.02
```

The Cartesian product of the supplied values is formed, but only combinations satisfying the strict Theorem-1 conditions (from the paper) are evaluated.

---

## 2. Distance-matrix version: `prune_deep_forest_distance.py`

### Input

One full, labeled, symmetric pairwise distance matrix in:

- TSV;
- CSV; or
- whitespace-delimited text.

Example TSV:

```text
taxon	A	B	C	D
A	0	0.10	0.40	0.50
B	0.10	0	0.42	0.51
C	0.40	0.42	0	0.12
D	0.50	0.51	0.12	0
```

The program checks that the matrix:

- contains at least 3 taxa;
- is square;
- has matching row and column labels;
- has unique, nonempty labels;
- contains only finite numeric values;
- has no negative distances;
- has a zero diagonal;
- is symmetric within numerical tolerance.

### Basic automatic-grid example

```bash
python prune_deep_forest_distance.py distance_matrix_8tips.tsv \
    --outdir prune_deep_output_8tips
```

For the provided 20-tip test matrix:

```bash
python prune_deep_forest_distance.py distance_matrix_20tips.tsv \
    --outdir prune_deep_output_distance_20tips
```

## Limit tau parameter refinement

```bash
python prune_deep_forest_distance.py distance_matrix_20tips.tsv \
    --outdir prune_deep_output_distance_20tips \
    --tau-refine-rounds 1
```

Disable refinement:

```bash
python prune_deep_forest_distance.py distance_matrix_20tips.tsv \
    --outdir prune_deep_output_20tips \
    --tau-refine-rounds 0
```

## Custom parameter grid

```bash
python prune_deep_forest_distance.py distance_matrix_20tips.tsv \
    --outdir prune_deep_output_20tips_custom \
    --M "[1.5,2.0]" \
    --m "[0.2,0.4,0.6]" \
    --tau "[0.005,0.02]"
```

As in the alignment-input version, combinations that do not satisfy the strict theorem conditions are discarded before Prune-Deep is evaluated.

---

## Automatic parameter search

The automatic search is the same in both scripts after a pairwise distance matrix has been obtained.

### 1. Neighbor Joining tree

NJ is inferred from the supplied/estimated pairwise distance matrix.

- In `prune_deep_forest.py`, this matrix is estimated from the alignment using JC69.
- In `prune_deep_forest_distance.py`, this matrix is supplied directly by the user.

NJ itself does not require JC69; it operates on the distance matrix it receives.

### 2. Tau search scale

Let the positive internal branch lengths of the NJ tree be

$$
b_1,\ldots,b_r.
$$

The automatic center is

$$
\tau_0=\frac{\operatorname{median}(b_1,\ldots,b_r)}{4}.
$$

The initial values are

$$
\tau=\tau_0\,[0.25,\;0.5,\;1,\;1.5,\;2].
$$

If the NJ tree has no positive internal branch, the implementation uses one quarter of the 10th percentile of positive pairwise distances as a fallback scale.

**Important:** this is an automatic-search heuristic. The theorem treats $\tau$ as a distance-distortion/error parameter; the NJ-based formula above is not itself a theorem-derived estimator of that error bound.

### 3. NJ chord-depth proxy and m

For every unique nontrivial split $A_e \mid B_e$ induced by an NJ internal edge, calculate

$$
c_e = \min_{i \in A_e,\; j \in B_e} \hat{d}(i,j).
$$

The NJ chord-depth proxy is

$$
\widehat{\Delta}_c^{NJ} = \max_e c_e.
$$

The automatic `m` grid is

$$
m = \widehat{\Delta}_c^{NJ} \times [0.5,\; 1.0,\; 1.5].
$$

This uses the same max-min chord-depth construction as the simulation benchmark, but replaces the unavailable true tree and true tree metric with the NJ topology and the observed/supplied distance matrix.

### 4. M values

For every $(m,\tau)$, define the theorem boundary

$$
M_{\text{base}}=2m+3\tau.
$$

The automatic candidates are

$$
M=M_{\text{base}}[1.05,\;1.25,\;1.50].
$$

All proposed triples are still checked against the full strict theorem conditions below.

---

## Theorem-1 parameter filtering

Only parameter combinations satisfying

$$
M>3\tau,\qquad m>3\tau,\qquad M>2m+3\tau
$$

are evaluated by Prune-Deep.

This filtering applies to both:

- automatic grids; and
- user-supplied custom grids.

Therefore a custom parameter combination can be supplied on the command line but still be discarded if it is outside the strict theorem-valid region.

---

## Adaptive lower-tau refinement

Adaptive tau refinement is used only with the automatic grid.

After the initial grid is evaluated:

1. select the current best split-compatible result;
2. check whether its $\tau$ is the smallest currently tested value;
3. if so, add

$$
\tau_{\text{new}}=\frac{\tau_{\min}}{2};
$$

4. generate theorem-valid $(M,m,\tau_{\text{new}})$ combinations;
5. evaluate them and reselect the best result;
6. continue until the winner is no longer at the lower tau boundary or the configured `--tau-refine-rounds` limit is reached.

Default:

```text
--tau-refine-rounds 3
```

Custom-grid runs do not automatically add smaller tau values.

---

## Reconstruction and selection

For each theorem-valid $(M,m,\tau)$ point, the program:

1. constructs the `m`-clustering graph using an edge when

   $$
   \hat d(i,j)<m;
   $$

2. treats connected components as forest components;
3. runs Mini Contractor only on edges of the `m`-clustering graph;
4. extends candidate local bipartitions by graph connectivity;
5. rejects ambiguous extensions;
6. canonicalizes and deduplicates recovered splits;
7. retains nontrivial splits, where both sides contain at least two taxa;
8. checks that all retained splits within a component are mutually compatible;
9. reconstructs a possibly unresolved Newick topology from the compatible split set.

For each result, define its resolution as the total number of recovered nontrivial splits. For a forest it is the sum of nontrivial splits. A single tree and a forest therefore compete using the same resolution measure.

The primary result is selected by:

1. maximum total nontrivial split count;
2. on an exact split-count tie, fewer components;
3. if still tied, earlier grid order.

All evaluated candidates already satisfy the theorem parameter filter.

---

## Output files

Each run writes a grid-search log and at least one Newick file.

### Primary result

The primary result is the globally most-resolved split-compatible candidate found by the search.

Default filename pattern:

```text
main_M{M}m{m}tau{tau}.nwk
```

If the result is a forest, the file contains one Newick record per component.

You can override the primary filename with:

```bash
--output my_result.nwk
```

### Secondary forest

If the primary result is a single tree and at least one compatible multi-component forest exists, the program also writes the most-resolved such forest:

```text
forest_less_preferred_M{M}m{m}tau{tau}.nwk
```

If the primary result is already a forest, no secondary forest file is written.

### Grid-search log

Default:

```text
grid_search.log
```

Override with:

```bash
--log my_search.log
```

The log includes:

- input information;
- NJ tree and positive internal branch lengths;
- automatic-grid scales or custom parameter values;
- theorem-filter counts;
- all evaluated grid points;
- split compatibility;
- number and sizes of forest components;
- total and component-wise nontrivial split counts;
- the selected result;
- optional secondary forest information;
- elapsed runtime.

For the alignment-input version, the log records:

```text
distance_model    JC69
```

For the distance-matrix version, the log records that the matrix was user supplied and that the distance model is unspecified by the program.

---

## Example test data

Two synthetic additive distance matrices are available for testing the distance-input implementation:

```text
distance_matrix_8tips.tsv
distance_matrix_20tips.tsv
```

Test them with:

```bash
python prune_deep_forest_distance.py distance_matrix_8tips.tsv \
    --outdir prune_deep_output_distance_8tips
```

and

```bash
python prune_deep_forest_distance.py distance_matrix_20tips.tsv \
    --outdir prune_deep_output_distance_20tips
```


---


© 2024. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.

This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
