import re
from io import StringIO
import random
import numpy as np
import networkx as nx
import dendropy
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix
np.random.seed(0)


def construct_clustering_graph(distance_matrix, M):
    num_leaves = distance_matrix.shape[0]
    G = nx.Graph()
    for i in range(num_leaves):
        G.add_node(i)
    for i in range(num_leaves):
        for j in range(i+1, num_leaves):
            if distance_matrix[i, j] < M:
                G.add_edge(i, j)

    return G

def extract_numbers(s):
    return re.findall(r'\d+', s)


def sort_and_freeze(pair):
    return tuple(sorted((frozenset(s) for s in pair), key=lambda s: (min(s) if s else float('inf'))))


#Th1 and 3.2
#m < 0.5(M - 3tau) and m >3tau; mini construction bipartition if edge length >= 4tau; d(u,v) < m
#Typically, M is much larger than tau. In that case, reconstruct a subforest of T with chord depth ≈1/2M which 
#includes all edges of length at least 4tau.
#all edges < tau will be contracted

# calculate phi in MC; this is a distance from node u to the intersection point of nodes (u, v, w)
def calculate_phi(u, v, w, distance_matrix):
    return 0.5 * (distance_matrix[u][v] + distance_matrix[u][w] - distance_matrix[v][w])

#Mini Contractor, works for pair of leaves, result depends on M, tau  
def mini_contractor(component, distance_matrix, leaves, M, tau, printout = 0):
    if printout == 0:
        u, v = leaves
        B = {w for w in component if max(distance_matrix[u][w], distance_matrix[v][w]) < M}
        #print("ball: ", B)
        if u not in B:
            B.add(u)  
        if v not in B:
            B.add(v)  
        S = B - {u}
        bipartitions = []
        x_minus_1 = u
        j = 0
        C = {0: [u]} #clusters
        phi = {w: calculate_phi(u, v, w, distance_matrix) for w in S}
        while S:
            x_0 = min(S, key=lambda w: phi[w])
            if j == 0: 
                phi[u] = 0
            if phi[x_0] - phi[x_minus_1] >= 2 * tau:
                bipartitions.append([B - S.copy(), S.copy()])
                C[j + 1] = [x_0]
                j += 1
            else:
                C[j].append(x_0)
            S.remove(x_0)
            x_minus_1 = x_0
    else:
        u, v = leaves
        B = {w for w in component if max(distance_matrix[u][w], distance_matrix[v][w]) < M}
        print("ball: ", B)
        if u not in B:
            B.add(u)  
        if v not in B:
            B.add(v)  
        print("**********************************")
        print("For a pair of leaves: ", u, " and ", " v ", v, " the ball is: ", B) 
        print("**********************************")
        S = B - {u}
        bipartitions = []
        x_minus_1 = u
        j = 0
        C = {0: [u]} #clusters
        phi = {w: calculate_phi(u, v, w, distance_matrix) for w in S}
        print("Phi from mini-contractor for each node w: ", phi)
        while S:
            x_0 = min(S, key=lambda w: phi[w])
            if j == 0: 
                phi[u] = 0
            print("Phi(x_0), Phi(x_minus_1), 2*tau: ", phi[x_0], phi[x_minus_1], 2*tau)
            if phi[x_0] - phi[x_minus_1] >= 2 * tau:
                print("long edge (bipartition) created between: ", x_0, " and ", x_minus_1, " diff: ", phi[x_0] - phi[x_minus_1])
                bipartitions.append([B - S.copy(), S.copy()])
                C[j + 1] = [x_0]
                j += 1
            else:
                print(x_0, " is close to ", x_minus_1, " no bipartition is created")
                C[j].append(x_0)
            print("Bipartitions and j: ", bipartitions, j)
            print("Clusters: ", C)
            # Update the set S (reducing) and the previous node
            S.remove(x_0)
            x_minus_1 = x_0
        print("Final Bipartitions and j: ", bipartitions, j)
        print("Final Clusters: ", C)
    return bipartitions, {k: list(v) for k, v in C.items()}


#Extender, 'graph' is the forest component 
def extender(graph, bipartitions, leaves, distance_matrix, printout=0):
    extended_bipartitions = []
    for psi_u, psi_v in bipartitions:
        #remove only cross edges 
        K = graph.copy()
        K.remove_edges_from(
            (a, b) for a in psi_u for b in psi_v if K.has_edge(a, b)
        )
        # extend by connectivity: every remaining leaf connects to exactly one side 
        ext_u, ext_v = set(), set()
        for comp in nx.connected_components(K):
            comp = set(comp)
            if comp & psi_u:
                ext_u |= comp
            elif comp & psi_v:
                ext_v |= comp
            else:
                pass
        extended_bipartitions.append((ext_u, ext_v))
    return extended_bipartitions


#pass bipartitions, get unique trees
def get_unique_trees(unique_sorted_tuples, Ntips, runs=10):
    unique_trees = set()
    taxon_namespace = dendropy.TaxonNamespace()

    for _ in range(runs):
        random.shuffle(unique_sorted_tuples)

        recovered_tree = reconstruct_unrooted_from_bipartitions(unique_sorted_tuples, Ntips)
        #print("tree: ", recovered_tree)
        #print('colless: ', treemeasure.colless_tree_imbalance(recovered_tree))
        recovered_tree_string = recovered_tree.as_string(schema="newick")
        recovered_newick_tree = re.search(r'(?=\()(.+;)', recovered_tree_string).group(1)

        is_unique = True
        for existing_tree_newick in unique_trees:
            existing_tree = dendropy.Tree.get(data=existing_tree_newick, schema="newick", taxon_namespace=taxon_namespace)
            new_tree = dendropy.Tree.get(data=recovered_newick_tree, schema="newick", taxon_namespace=taxon_namespace)
            rf_distance = dendropy.calculate.treecompare.symmetric_difference(existing_tree, new_tree)
            if rf_distance == 0:
                is_unique = False
                break

        if is_unique:
            unique_trees.add(recovered_newick_tree)

    return unique_trees


def get_internal_edge_lengths(tree):
    internal_edge_lengths = []
    for node in tree.traverse():
        if not node.is_leaf() and node.dist is not None:
        #if node.dist is not None:
            internal_edge_lengths.append(node.dist)
    return internal_edge_lengths


def create_mappings(unique_sorted_tuples):
    unique_numbers = set()
    for A, B in unique_sorted_tuples:
        unique_numbers.update(A, B)
    sorted_numbers = sorted(unique_numbers)
    # mapping from original to continuous range
    mapping = {number: i for i, number in enumerate(sorted_numbers)}
    # reverse mapping from continuous range back to original
    reverse_mapping = {i: number for number, i in mapping.items()}
    return mapping, reverse_mapping

def reconstruct_unrooted_from_bipartitions(unique_sorted_tuples, num_taxa):
    #print('unique_sorted_tuples: ', unique_sorted_tuples, 'num_taxa ', num_taxa)
    mapping, reverse_mapping = create_mappings(unique_sorted_tuples)
    number_strings = [str(i) for i in range(num_taxa)]
    labels_names = dendropy.TaxonNamespace(number_strings, label="taxa1")
    split_matrix = np.zeros((len(unique_sorted_tuples), num_taxa))
    #print("mapping, reverse_mapping")
    #print(mapping, reverse_mapping)

    for i, (A, B) in enumerate(unique_sorted_tuples):
        for a in A:
            mapped_a = mapping[a]
            split_matrix[i, mapped_a] = 1
    #print(split_matrix.astype(int))
    bitmask_splits = []
    for bipartition in split_matrix:
        bitmask = 0
        for index, value in enumerate(bipartition):
            if value == 1:
                bitmask |= 1 << index
        bitmask_splits.append(bitmask)
    #print(bitmask_splits)
    number_strings = [str(i) for i in range(num_taxa)]
    #tns1 = dendropy.TaxonNamespace(number_strings, label="taxa1")
    my_unrooted_tree = dendropy.Tree.from_split_bitmasks(bitmask_splits, taxon_namespace = labels_names)
    #print(my_unrooted_tree.as_string(schema="newick"))
    for taxon in my_unrooted_tree.taxon_namespace:
        original_number = reverse_mapping[int(taxon.label)]
        taxon.label = str(original_number)
    #print(my_unrooted_tree.as_string(schema="newick"))

    return my_unrooted_tree

def find_identical_species(distance_matrix):
    identical_species = False
    for i in range(distance_matrix.shape[0]):
        for j in range(i + 1, distance_matrix.shape[1]):
            if distance_matrix.iloc[i, j] == 0:
                print(f"Species {distance_matrix.index[i]} and {distance_matrix.columns[j]} are identical.")
                identical_species = True
    if not identical_species:
        print("No identical species found.")
        return 1
    else:
        return 0
        

def jc69_distance(sequence1, sequence2):
    n = len(sequence1)
    diffs = sum(1 for x, y in zip(sequence1, sequence2) if x != y)
    p = diffs / n
    if p < 0.75:
        distance = -3/4 * np.log(1 - 4/3 * p)
    else:
        distance = float('inf')  
    return distance

def calculate_jc69_distance_matrix(msa):
    num_seqs = len(msa)
    labels = [record.id for record in msa]
    matrix = []
    for i in range(num_seqs):
        row = []
        for j in range(i + 1):  
            if i == j:
                row.append(0.0)  
            else:
                dist = jc69_distance(str(msa[j].seq), str(msa[i].seq))
                row.append(dist)
        matrix.append(row)
    return DistanceMatrix(names=labels, matrix=matrix)




def calculate_chord_depth(newick_string, distance_matrix):
    tree = dendropy.Tree.get(data=newick_string, schema="newick")
    max_min_path_length = 0
    def collect_leaf_labels(node):
        return [leaf.taxon.label for leaf in node.leaf_iter()]
    for edge in tree.postorder_edge_iter():
        if not edge.is_terminal():
            min_path_length = float('inf')            
            child_node = edge.head_node
            parent_node = edge.tail_node
            if parent_node and child_node:
                parent_node.remove_child(child_node)
                component1 = set(collect_leaf_labels(child_node))
                component2 = set(collect_leaf_labels(tree.seed_node)) - component1  # seed_node is the root node of the tree
                parent_node.add_child(child_node)
                for taxon1 in component1:
                    for taxon2 in component2:
                        try:
                            path_length = distance_matrix[taxon1][taxon2]
                            if path_length < min_path_length:
                                min_path_length = path_length
                        except KeyError:
                            print(f"Missing distance data for taxa: {taxon1}, {taxon2}")

            if min_path_length < float('inf') and min_path_length > max_min_path_length:
                max_min_path_length = min_path_length
    
    return max_min_path_length




def prune_tree_and_remove_branch_lengths(tree, species_list):
    tree.prune(species_list, preserve_branch_length=False)
    return tree.write(format=9)



def adjust_taxon_labels(newick):
    def replace_label(match):
        num = int(match.group(1))  
        return str(num - 1)  
    adjusted_newick = re.sub(r'T(\d+)', replace_label, newick)
    return adjusted_newick


def adjust_labels(labels):
    new_labels = []
    for label in labels:
        if label.startswith('T') and label[1:].isdigit():
            new_label = str(int(label[1:]) - 1)
        else:
            new_label = label  # keep the label unchanged if it doesn't match expected format
        new_labels.append(new_label)
    return new_labels



def count_internal_nodes(newick_string):
    tree = Phylo.read(StringIO(newick_string), 'newick')
    internal_nodes = [node for node in tree.get_nonterminals()]
    return len(internal_nodes)


def calculate_internal_node_ratio(int_nodes, component_sizes, num_components):
    if num_components > 1:
        internal_node_ratios = [int_node / component_size for int_node, component_size in zip(int_nodes, component_sizes)]
    else:
        internal_node_ratios = [int_nodes[0] / component_sizes[0]]
    return internal_node_ratios



# formatting helper for tabular TXT (TSV) output
def _fmt(value):
    if value is None:
        return "NA"
    if isinstance(value, (list, tuple)):
        return "[" + ",".join(str(x) for x in value) + "]"
    try:
        import numpy as _np
        if isinstance(value, _np.generic):
            return str(value.item())
        if isinstance(value, _np.ndarray):
            return "[" + ",".join(str(x) for x in value.tolist()) + "]"
    except Exception:
        pass
    return str(value)

def _fmt_list(x):
    try:
        arr = np.asarray(x, dtype=float).tolist()
        return "[" + ",".join(f"{v:.6g}" for v in arr) + "]"
    except Exception:
        return str(x)