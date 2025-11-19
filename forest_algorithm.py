#also requires IQTree2 executable or already simulated alignment file (see main())
import re
import io
import os, subprocess
import argparse
import numpy as np
import networkx as nx
import pandas as pd
import dendropy
from Bio import Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from ete3 import Tree
from itertools import product
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
np.random.seed(0)

#prune_deep_functions.py contains helper functions
from prune_deep_helper_functions import (
    construct_clustering_graph, extract_numbers, sort_and_freeze, 
    mini_contractor, extender, get_unique_trees, 
    get_internal_edge_lengths, 
    find_identical_species, calculate_jc69_distance_matrix, 
    calculate_chord_depth, prune_tree_and_remove_branch_lengths, adjust_taxon_labels, adjust_labels, 
    count_internal_nodes, calculate_internal_node_ratio, _fmt, _fmt_list
)

#Runs splits_compatible.R file that takes 2 trees and tells whether all the splits in tree2 are also in tree1 (pruned to
#the tips in tree2). Also outputs #splits in tree2
def split_R(tree1_newick, tree2_newick):
    result = subprocess.run(
        ["/usr/projects/t23_ad6800/condaenvs/phylo_env/bin/Rscript", "splits_compatible.R", tree1_newick, tree2_newick],
        capture_output=True, text=True
    )
    output_lines = result.stdout.strip().splitlines()
    
    if output_lines and output_lines[0].startswith('TRUE'):
        is_compatible = True
    elif output_lines and output_lines[0].startswith('FALSE'):
        is_compatible = False
    else:
        raise ValueError("Unexpected output from R script:", output_lines)
    
    split_info = output_lines[0].split()
    if len(split_info) > 1:
        split_count = int(split_info[1])
    else:
        raise ValueError("Expected a number after TRUE/FALSE, but none was found:", output_lines)
    
    return is_compatible, split_count

#Same function as split_R but also outputs false splits if any
def false_splits_R(tree1_newick, tree2_newick):
    result = subprocess.run(
        ["/usr/projects/t23_ad6800/condaenvs/phylo_env/bin/Rscript", "false_splits.R", tree1_newick, tree2_newick],
        capture_output=True, text=True
    )
    output_lines = result.stdout.strip().splitlines()

    is_compatible = None
    split_count = None
    false_splits = []

    for i, line in enumerate(output_lines):
        if line.startswith('$istrue'):
            is_compatible = output_lines[i + 1].strip() == '[1] "TRUE"'
        elif line.startswith('$numsplits'):
            split_count = int(output_lines[i + 1].strip().split()[1].replace('"', ''))
        elif line.startswith('$falsesplits'):
            false_splits = int(output_lines[i + 1].split(']', 1)[-1].strip().strip('"').split()[0])

    if is_compatible is None:
        raise ValueError("Could not determine compatibility from output:", output_lines)
    if split_count is None:
        raise ValueError("Could not determine number of splits from output:", output_lines)
    return is_compatible, split_count, false_splits




# Algorithm function (kept same, returns list of dict results) 
# input: 1/2*shortest branch length in the true tree, true tree newick string, NJ newick string, true distance matrix (pandas df), estimated distance matrix (pandas df), values of ms, taus, Ms, number of tips, 
# threshold for min acceptable #nodes in forest component, printout variable to hide/show full output in terminal
def get_forest_split_check(min_len_true, true_newick_string, nj_tree, true_a, a, ms, taus, Ms, Ntips, Thr, printout=0):
    collected_data = []
    t = Tree(true_newick_string)
    njt = Tree(nj_tree)
    chord_d = calculate_chord_depth(true_newick_string, true_a)
    grid_size = 0
    if find_identical_species(a) == 1:    
        for m in ms:
            for tau in taus:
                for M in Ms:
                    if M > 2*m + 3*tau and m > 3*tau: #Theorem 1 condition
                        true_g = construct_clustering_graph(np.array(a), m)
                        num_components = nx.number_connected_components(true_g)
                        components = list(nx.connected_components(true_g)) #split graph in components, depends on m
                        if all(len(component) > Thr for component in components):
                            all_components_single_tree = True
                            component_trees = {}
                            component_sizes = []
                            for i, component in enumerate(components): #for each component run Forest algorithm
                                component_subgraph = true_g.subgraph(component)
                                component_size = len(component_subgraph.nodes())
                                component_sizes.append(component_size)
                                leaves_pairs = [(u, v) for u in component_subgraph.nodes() for v in component_subgraph.nodes() if u < v]
                                final_bipartitions = []

                                for leaves in leaves_pairs:
                                    mini_bipartitions, clusters = mini_contractor(component_subgraph, np.array(a), leaves, M, tau, printout=printout)
                                    extended_bipartition = extender(component_subgraph, mini_bipartitions, leaves, np.array(a))
                                    final_bipartitions.append(extended_bipartition)
                                    
                                flattened_bipartitions = [tuple(frozenset(part) for part in bipartition) for sublist in final_bipartitions for bipartition in sublist]
                                unique_bipartitions = list(set(flattened_bipartitions))
                                unique_bipartitions_as_sets = [tuple(set(part) for part in bipartition) for bipartition in unique_bipartitions]
                                sorted_frozenset_data = [sort_and_freeze(pair) for pair in unique_bipartitions_as_sets]
                                unique_sorted_frozenset_data = set(sorted_frozenset_data)
                                unique_sorted_tuples = [tuple(map(set, pair)) for pair in unique_sorted_frozenset_data]
                                #tree popping to get trees from bipartitions
                                if len(unique_sorted_tuples) > 0:
                                    unique_trees = get_unique_trees(unique_sorted_tuples, len(component_subgraph.nodes()), 10)
                                    if len(unique_trees) == 1:
                                        component_trees[i] = unique_trees
                                    else:
                                        all_components_single_tree = False
                                        break
                                else:
                                    all_components_single_tree = False
                                    break
                            #only care if we have 1 unique tree resulting from bipartitions
                            if all_components_single_tree: 
                                induced_rfs_list = []
                                tree_parts = []
                                induced_rfs_nj_list = []
                                nj_parts = []
                                int_nodes = []
                                res_nodes = []
                                is_compatible_list = []
                                split_count_list = []
                                num_false_splits_list = []
                                grid_size += 1
                                #for each forest component prune true and NJ trees to taxa in the component and 
                                #calculate RF distance, whether splits are compatible; and number of false splits if any
                                for component_idx, trees in component_trees.items():
                                    for sub_newick in trees:
                                        taxon_namespace = dendropy.TaxonNamespace()
                                        forest_part_tree = dendropy.Tree.get(data=sub_newick, schema="newick", taxon_namespace=taxon_namespace)
                                        tree_parts.append(sub_newick)

                                        # species list for pruning
                                        species_list = extract_numbers(sub_newick)

                                        # true induced part
                                        true_tree_cp = t.copy()
                                        true_tree_cp.prune(species_list, preserve_branch_length=False)
                                        true_part_tree = dendropy.Tree.get(data=true_tree_cp.write(), schema="newick", taxon_namespace=taxon_namespace)

                                        # RF: forest vs true (induced)
                                        induced_rf_val = dendropy.calculate.treecompare.symmetric_difference(true_part_tree, forest_part_tree)
                                        induced_rfs_list.append(induced_rf_val)
                                        
                                        
                                        # true vs NJ induced part
                                        int_nodes.append(count_internal_nodes(sub_newick))
                                        nj_cp = njt.copy()
                                        nj_induced_newick = prune_tree_and_remove_branch_lengths(nj_cp, species_list)
                                        nj_parts.append(nj_induced_newick)
                                        nj_part_tree = dendropy.Tree.get(data=nj_induced_newick, schema="newick", taxon_namespace=taxon_namespace)
                                        induced_rf_nj_val = dendropy.calculate.treecompare.symmetric_difference(true_part_tree, nj_part_tree)
                                        induced_rfs_nj_list.append(induced_rf_nj_val)  

                                        res_nodes.append(sub_newick.count('('))
                                        
                                        # Split compatibilities
                                        is_compat, split_cnt = split_R(true_newick_string, sub_newick)
                                        is_compatible_list.append(is_compat)
                                        split_count_list.append(split_cnt)
                                        
                                        # False splits 
                                        _, _, num_false_splits_val = false_splits_R(true_newick_string, sub_newick) 
                                        num_false_splits_list.append(num_false_splits_val)
                                #dictionary of the result run        
                                component_info = {
                                    'M': M,
                                    'm': m,
                                    'tau': tau,
                                    'grid_size': grid_size, #size of the grid if any
                                    'min_len_true': min_len_true, #shortest branch length in the true tree
                                    'max_dhat': a.max().max(), #max value in the estimated distance matrix
                                    'max_d': true_a.max().max(), #max value in the true distance matrix
                                    'chord_true_depth': chord_d, #chord depth of the true tree
                                    'd_less_Mplustau': 1 if true_a.max().max() < M + tau else 0, #distortion condition check (see paper)
                                    'd_dhat_less_Mplustau': 1 if a.max().max() < M + tau else 0, #distortion condition check
                                    'd_minus_dhat_less_tau': 1 if np.abs(true_a - a).max().max() < tau else 0, #distortion condition check
                                    'd_diff': np.abs(true_a - a).max().max(), 
                                    '%_diff_greater_tau': np.sum(np.sum(np.abs(true_a - a) > tau))/(Ntips*Ntips - Ntips)*100, #number of pairs for which |d-d^|> tau (violation of the distortion)
                                    'num_components': num_components, #number of forest components
                                    'component_sizes': component_sizes if num_components > 1 else component_sizes[0], #forest component sizes
                                    'induced_rfs': induced_rfs_list if num_components > 1 else [induced_rfs_list[0]], #RF(true, component) for each component
                                    'resolved_nodes': np.sum(res_nodes) if num_components > 1 else [res_nodes[0]], #how many nodes resolved
                                    'is_compatible': is_compatible_list if num_components > 1 else [is_compatible_list[0]], #whether splits in component are compatible
                                    'split_count': split_count_list if num_components > 1 else [split_count_list[0]],#count #splits for each component
                                    'sum_split_count': np.sum(split_count_list) if num_components > 1 else [split_count_list[0]], #sum above (#splits)
                                    'forest_num_falsesplit_count': num_false_splits_list, #number of false splits if any
                                    'forest_sum_falsesplit_count': float(np.sum(num_false_splits_list)) if num_components > 1 else [float(num_false_splits_list[0])], #sum above across components
                                    'number_internal_nodes': int_nodes if num_components > 1 else [int_nodes[0]], #number of internal (not terminal) nodes in forest components
                                    'ratio_internal_nodes': np.round(np.mean(calculate_internal_node_ratio(int_nodes, component_sizes, num_components)), 2), #weighted #internal nodes by component sizes
                                    'forest': tree_parts if num_components > 1 else tree_parts[0], #forest tree components
                                    'induced_rfs_nj': induced_rfs_nj_list if num_components > 1 else [induced_rfs_nj_list[0]], #RF(true, NJ)
                                    'nj_parts': nj_parts if num_components > 1 else nj_parts[0], #NJ tree components (pruned to match forest taxa)
                                    'NJ_#internal_nodes': count_internal_nodes(nj_tree) #number of internal (not terminal) nodes in NJ components
                                }
                                if printout:
                                    print(component_info)
                                collected_data.append(component_info)
                            else:
                                component_info = {
                                    'M': M,
                                    'm': m,
                                    'tau': tau,
                                    'max_dhat': a.max().max(),
                                    'max_d': true_a.max().max(),
                                    'chord_true_depth': chord_d,
                                    'd_less_Mplustau': 1 if true_a.max().max() < M + tau else 0,
                                    'd_dhat_less_Mplustau': 1 if a.max().max() < M + tau else 0,
                                    'd_minus_dhat_less_tau': 1 if np.abs(true_a - a).max().max() < tau else 0,
                                    'num_components': num_components,
                                    'component_sizes': component_sizes if num_components > 1 else component_sizes[0],
                                    'induced_rfs': None,
                                    'forest_polydist': None,
                                    'forest': None,
                                    'induced_rfs_nj': None,
                                    'nj_parts': None
                                }
                                collected_data.append(component_info)
    else:
        print("Identical species found.") #happens if sequence length k is too low
    return collected_data


# One-process worker to evaluate a single (m, M, tau) grid combination
def _run_one_combo(rep_id,
                   m, M, tau,
                   min_len_true,
                   true_newick_string,
                   nj_tree_string,
                   true_distance_matrix,
                   est_dist_mat,
                   ntips, Thr):
    try:
        res_list = get_forest_split_check(
            min_len_true=min_len_true,
            true_newick_string=true_newick_string,
            nj_tree=nj_tree_string,
            true_a=true_distance_matrix,
            a=est_dist_mat,
            ms=[m], taus=[tau], Ms=[M],
            Ntips=ntips, Thr=Thr, printout=0
        )
        if not res_list:
            return {
                "replicate": rep_id, "m": m, "M": M, "tau": tau,
                "num_components": "NA",
                "component_sizes": "NA",
                "is_compatible": "NA",
                "split_count": "NA",
                "sum_split_count": "NA",
                "forest_num_falsesplit_count": "NA",
                "forest_sum_falsesplit_count": "NA",
                "pct_diff_greater_tau": "NA",
            }
        c = res_list[0]
        row = {
            "replicate": rep_id,
            "m": m, "M": M, "tau": tau,
            "num_components": _fmt(c.get("num_components")),
            "component_sizes": _fmt(c.get("component_sizes")),
            "is_compatible": _fmt(c.get("is_compatible")),
            "split_count": _fmt(c.get("split_count")),
            "sum_split_count": _fmt(c.get("sum_split_count")),
            "forest_num_falsesplit_count": _fmt(c.get("forest_num_falsesplit_count")),
            "forest_sum_falsesplit_count": _fmt(c.get("forest_sum_falsesplit_count")),
            "pct_diff_greater_tau": _fmt(c.get("%_diff_greater_tau")),
        }
        return row
    except Exception as e:
        # return an error row (so the output file table stays aligned)
        return {
            "replicate": rep_id, "m": m, "M": M, "tau": tau,
            "num_components": f"ERROR:{e}",
            "component_sizes": "NA",
            "is_compatible": "NA",
            "split_count": "NA",
            "sum_split_count": "NA",
            "forest_num_falsesplit_count": "NA",
            "forest_sum_falsesplit_count": "NA",
            "pct_diff_greater_tau": "NA",
        }

# Parse k and n from command line; allow one or more integers for k(sequence length) and n (#tips)
def parse_args():
    p = argparse.ArgumentParser()
    
    p.add_argument("-k", "--k", "--k-values",
                   dest="k_values", type=int, nargs="+", default=[512],
                   help="One or more k values (e.g., --k 512 or --k 512 1024)")
    p.add_argument("-n", "--n", "--Ntips",
                   dest="N_tips", type=int, nargs="+", default=[128],
                   help="One or more k values (e.g., --n 100 or --n 128 256)")  
    return p.parse_args()




# Main 
def main():
    print("Starting script - Reading existing files")

    # ensure working dir is script dir
    os.chdir(os.path.dirname(os.path.abspath(__file__)))

    #user feeds k and n 
    np.random.seed(2024)
    args = parse_args()
    k_values = args.k_values
    Ntips = args.N_tips 
    ntips = Ntips[0]
    Thr = 2
    num_aln_reps = 2 #how many alignments to simulate via IQTree2

    forest_results_128 = {}
    fin_res = []

    # read existing true tree
    true_tree_file = f"true_tree_{ntips}.tre"
    print(f"Reading tree from: {true_tree_file}")
    with open(true_tree_file, 'r') as f:
        true_tree_content = f.read().strip()

    taxon_namespace = dendropy.TaxonNamespace()
    true_tree_string = adjust_taxon_labels(true_tree_content)
    print(f"True tree string: {true_tree_string}")

    # internal edges
    internal_edges = get_internal_edge_lengths(Tree(true_tree_string))
    print(f"internal_edges: {np.sort(internal_edges)}")
    print(f"0.5*internal_edges: {0.5*np.sort(internal_edges)}")

    true_tree = dendropy.Tree.get(data=true_tree_string, schema="newick", taxon_namespace=taxon_namespace)

    # distance matrix from true tree
    pdm = dendropy.calculate.treemeasure.PatristicDistanceMatrix(true_tree)
    labels = [taxon.label for taxon in pdm.taxon_namespace]
    matrix_size = len(pdm.taxon_namespace)
    distance_matrix = [[0 for _ in range(matrix_size)] for _ in range(matrix_size)]
    for i, taxon1 in enumerate(pdm.taxon_namespace):
        for j, taxon2 in enumerate(pdm.taxon_namespace):
            distance_matrix[i][j] = pdm(taxon1, taxon2)
    true_distance_matrix = pd.DataFrame(distance_matrix, index=labels, columns=labels)
    true_distance_matrix.index = adjust_labels(true_distance_matrix.index)
    true_distance_matrix.columns = adjust_labels(true_distance_matrix.columns)
    true_distance_matrix = true_distance_matrix.sort_index().sort_index(axis=1)
    true_distance_matrix.index = true_distance_matrix.index.astype(int)
    true_distance_matrix = true_distance_matrix.sort_index()
    true_distance_matrix.index = true_distance_matrix.index.astype(str)
    true_distance_matrix.columns = true_distance_matrix.columns.astype(int)
    true_distance_matrix = true_distance_matrix.sort_index(axis=1)
    true_distance_matrix.columns = true_distance_matrix.columns.astype(str)

    chord_true_depth_value = calculate_chord_depth(true_tree_string, true_distance_matrix)
    half_internal_edges = 0.5 * np.sort(internal_edges).astype(float)

    # paths to working directory and IQTree executable
    tipdir = "./distphylo"
    iqtree2_path = "./distphylo/iqtree2" 
    treefile = f"true_tree_{ntips}.tre"

    # grid choices for (m, M, tau) for this run; USER can change
    Ms = [round(2.5 - 0.1*i, 1) for i in range(int((2.5-0.5)/0.1) + 1)] #[2.5, 2.2, 2, 1.75, 1.5, 1.25, 1, 0.75]
    ms = [0.8, 0.7, 0.65, 0.6, 0.55, 0.5, 0.45, 0.4]
    taus = np.linspace(0.01, 0.10, num=19).tolist() #[0.12, 0.11, 0.1, 0.095, 0.09, 0.085, 0.08, 0.075, 0.07, 0.065, 0.06, 0.055, 0.05, 0.045, 0.04]
    # work per sequence length k
    for k in k_values:
        print(f"seq_length: {k}")
        alnfile_prefix = os.path.join(tipdir, f"aln_{ntips}_{k}")
        print(alnfile_prefix)

        # simulate alignments via IQ-TREE2 (writes aln_k_i.fa to tipdir) or comment the following 3 lines to just read them from tipdir
        # or install IQTree with conda    
        cmd = [iqtree2_path, "--alisim", alnfile_prefix, "-t", treefile, "-m", "JC",
               "--length", str(k), "--num-alignments", str(num_aln_reps), "--out-format", "fasta"]
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        #write down replicates # to work with, if num_aln_reps = 2 then (1,2) - two replicates
        for i_replicate in (1, 2):
            print(f"replicate: {i_replicate}")
            alignment_path = os.path.join(tipdir, f"aln_{ntips}_{k}_{i_replicate}.fa")
            print(f"Reading alignment from: {alignment_path}")

            # build JC69 distance matrix from MSA
            alignments = SeqIO.parse(alignment_path, 'fasta')
            msa = MultipleSeqAlignment(alignments)
            est_dist_mat_orig = calculate_jc69_distance_matrix(msa)

            matrix_data = est_dist_mat_orig.matrix
            labels = est_dist_mat_orig.names
            full_matrix = [[0] * len(labels) for _ in range(len(labels))]
            for i in range(len(labels)):
                for j in range(i + 1):
                    full_matrix[i][j] = matrix_data[i][j]
                    full_matrix[j][i] = matrix_data[i][j]
            est_dist_mat = pd.DataFrame(full_matrix, index=labels, columns=labels) #estimated distance matrix d-hat

            # NJ tree from estimated distances
            constructor = DistanceTreeConstructor()
            tree = constructor.nj(est_dist_mat_orig)
            nj_tree_buf = io.StringIO()
            Phylo.write(tree, nj_tree_buf, "newick")
            nj_tree_string = re.sub(r'Inner\d+:', '', nj_tree_buf.getvalue())
            nj_tree_string = adjust_taxon_labels(nj_tree_string)
            nj_tree = dendropy.Tree.get(data=nj_tree_string.rstrip('\n'), schema="newick", taxon_namespace=taxon_namespace)
            rf_distance = dendropy.calculate.treecompare.symmetric_difference(true_tree, nj_tree)
            print(f'rf_distance: {rf_distance}')

            print(f"0.5*internal_edges: {0.5*np.sort(internal_edges)}")


            # Parallel grid over (m, M, tau) 
            # filter grid combos by Theorem 1 constraints
            combos = [(m, M, tau) for (m, M, tau) in product(ms, Ms, taus)
                      if (M > 2*m + 3*tau) and (m > 3*tau)]

            # 1/2*shortest branch length in the true tree
            min_len_true = float(0.5 * np.sort(internal_edges)[0])

            # how many workers to spawn
            workers_env = os.environ.get("SLURM_CPUS_ON_NODE")
            if workers_env and workers_env.isdigit():
                max_workers = max(1, int(workers_env))
            else:
                max_workers = max(1, mp.cpu_count() - 1)

            rep_rows = []
            with ProcessPoolExecutor(max_workers=max_workers) as ex:
                futs = []
                for m_val, M_val, tau_val in combos:
                    fut = ex.submit(
                        _run_one_combo,
                        i_replicate,  # replicate id
                        m_val, M_val, tau_val,
                        min_len_true,
                        true_tree_string,
                        nj_tree_string,
                        true_distance_matrix,  
                        est_dist_mat,         
                        ntips, Thr
                    )
                    futs.append(fut)

                for fut in as_completed(futs):
                    rep_rows.append(fut.result())

            # store raw results per replicate
            key = f"forest_results_{i_replicate}"
            forest_results_128[key] = rep_rows  
            # a global list for final TXT
            fin_res.extend(rep_rows)

            print(f"Finished replicate {i_replicate}: {len(rep_rows)} grid rows")
            print("-------------------------------------------------")

            # write one TXT table for all replicates 
            # out_txt = f"grid_summary_ntips{ntips}_{i_replicate}_k{'_'.join(map(str, k_values))}.txt"
            # header = [
            #     "replicate", "m", "M", "tau",
            #     "num_components", "component_sizes",
            #     "is_compatible", "split_count", "sum_split_count",
            #     "forest_num_falsesplit_count", "forest_sum_falsesplit_count",
            #     "pct_diff_greater_tau"
            # ]
            # with open(out_txt, "w") as fh:
            #     fh.write("\t".join(header) + "\n")
            #     for row in fin_res:
            #         fh.write("\t".join(_fmt(row[h]) for h in header) + "\n")

            # print(f"Wrote table: {out_txt}")


            # filter (no NAs/ERRORs), sort, and write per-replicate TSV 
            header_cols = [
                "replicate", "m", "M", "tau",
                "num_components", "component_sizes",
                "is_compatible", "split_count", "sum_split_count",
                "forest_num_falsesplit_count", "forest_sum_falsesplit_count",
                "pct_diff_greater_tau"
            ]
            df = pd.DataFrame(rep_rows)

            # all expected columns exist (fill missing with NaN)
            for col in header_cols:
                if col not in df.columns:
                    df[col] = np.nan
            df = df.replace(to_replace=r'^NA$', value=np.nan, regex=True)
            df = df.replace(to_replace=r'^ERROR:.*$', value=np.nan, regex=True)
            for col in ["m", "tau", "M"]:
                df[col] = pd.to_numeric(df[col], errors="coerce")

            # drop rows with ANY NaN in the grid columns for readability 
            df_clean = df[header_cols].dropna(how="any")

            # sort by m, then tau, then M (all descending) for convenience
            df_sorted = df_clean.sort_values(by=["m", "tau", "M"],
                                            ascending=[False, False, False],
                                            kind="mergesort")

            # output path (sorted)
            out_tsv = f"grid_summary_ntips{ntips}_{i_replicate}_k{k}_sorted.tsv"
            with open(out_tsv, "w") as fh:
                fh.write(f"chord_true_depth\t{chord_true_depth_value}\n")
                fh.write(f"0.5*internal_edges\t{_fmt_list(half_internal_edges)}\n")
                fh.write(f"seq_length\t{k}\n")
                fh.write(f"ntips\t{ntips}\n")
                fh.write(f"k\t{k}\n")
                fh.write(f"i_replicate\t{i_replicate}\n")
                fh.write(f"alignment_path\t{alignment_path}\n")
                fh.write(f"rf_distance\t{rf_distance}\n")
                fh.write("\n")
                df_sorted.to_csv(fh, sep="\t", index=False)
            print(f"Wrote sorted (no-NA) TSV with metadata: {out_tsv}")

if __name__ == "__main__":
    main()