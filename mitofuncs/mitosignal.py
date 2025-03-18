import os
import Bio
import random
import itertools
import numpy as np
import pandas as pd
import seaborn as sn
import matplotlib.pyplot as plt
from pathlib import Path
from Bio import Entrez,SeqIO
import time
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, to_tree, leaves_list
import ete3
from ete3 import Tree
from collections import defaultdict
import networkx as nx

def damerau_levenshtein_distance(genelist_1, genelist_2):
    """
    Link to wikipedia page: https://en.wikipedia.org/wiki/Damerau-Levenshtein_distance
    """
    # Define alphabet size based on unique elements in genelist_1 and genelist_2
    sigma = set(genelist_1) | set(genelist_2)
    da = defaultdict(int)  # Dictionary to store last seen positions
    
    len_a, len_b = len(genelist_1), len(genelist_2)
    maxdist = len_a + len_b
    
    # Initialize distance matrix
    d = [[maxdist] * (len_b + 2) for _ in range(len_a + 2)]
    
    for i in range(len_a + 1):
        d[i + 1][0] = maxdist
        d[i + 1][1] = i
    for j in range(len_b + 1):
        d[0][j + 1] = maxdist
        d[1][j + 1] = j
    
    for i in range(1, len_a + 1):
        db = 0
        for j in range(1, len_b + 1):
            k = da[genelist_2[j - 1]]
            l = db
            
            if genelist_1[i - 1] == genelist_2[j - 1]:
                cost = 0
                db = j
            else:
                cost = 1
            
            d[i + 1][j + 1] = min(
                d[i][j] + cost,      # Substitution
                d[i + 1][j] + 1,     # Insertion
                d[i][j + 1] + 1,     # Deletion
                d[k][l] + (i - k - 1) + 1 + (j - l - 1)  # Transposition
            )
            
        da[genelist_1[i - 1]] = i
    
    return d[len_a + 1][len_b + 1]

def double_cut_and_join_distance(genelist_1, genelist_2):
    """
    Link to wikipedia page: https://en.wikipedia.org/wiki/Double_Cut_and_Join_Model
    """
    # Create adjacency graphs from gene orders
    def build_adjacency_graph(genelist):
        adj_graph = {}  # Dictionary representing adjacency pairs
        for i in range(len(genelist) - 1):
            adj_graph[genelist[i]] = genelist[i + 1]
            adj_graph[genelist[i + 1]] = genelist[i]
        return adj_graph
    
    adj_graph_1 = build_adjacency_graph(genelist_1)
    adj_graph_2 = build_adjacency_graph(genelist_2)
    
    # Construct a breakpoint graph as an undirected graph
    breakpoint_graph = nx.Graph()
    for gene in adj_graph_1:
        if gene in adj_graph_2:
            breakpoint_graph.add_edge(adj_graph_1[gene], adj_graph_2[gene])
    
    num_cycles = nx.number_connected_components(breakpoint_graph)
    num_genes = len(genelist_1)
    
    # DCJ-distance formula: d = num_genes - num_cycles
    return num_genes - num_cycles

def gene_order_distance(genelist_1, genelist_2, metric_name):
    if metric_name=="damerau_levenshtein":
        return damerau_levenshtein_distance(genelist_1, genelist_2)
    elif metric_name=="double_cut_and_join":
        return double_cut_and_join_distance(genelist_1, genelist_2)
    
def arrange_in_order_based_on_distances(gene_orders_distances_dataframe):
    # Perform hierarchical clustering
    Z = linkage(gene_orders_distances_dataframe, method='average')

    # Get the order of the leaves
    ordered_indices = leaves_list(Z)

    # Get the species names in the new order
    evolutionary_order = gene_orders_distances_dataframe.index[ordered_indices].tolist()

    return evolutionary_order