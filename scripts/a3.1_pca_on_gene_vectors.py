"""
THIS FILE IS MEANT TO GENERATE PLOTS (VERY ROUGH IDEAS) FOR THE 
VARIOUS POSISBLE TRENDS IN THE GENE COUNTS AND GENE ORDERS
"""

from mitofuncs.mitoevo import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np
from sklearn.decomposition import PCA

all_species_names = list_all_species_names_from_file_path()
all_species_genomes = {species: parse_mitochondrial_genome(f"data/genbank_files/{species}/{species}_mitochondrion.gb") for species in all_species_names}
all_species_gene_dataframes = {species: pd.read_csv(f"data/genbank_files/{species}/{species}_cleaned_gene_data.tsv", delimiter="\t").squeeze() for species in all_species_names}

gene_trends_folder_path = Path("results/gene_trends")
if not gene_trends_folder_path.exists():
    gene_trends_folder_path.mkdir(parents=True)

"""
Since you have gene orders for multiple chromosomes and varying gene counts, you need a way to embed these gene orders into a fixed-dimensional space. Here are a few approaches you can take:

1. k-mer Based Frequency Vectors
Construct k-mers (subsequences of length k) from each chromosome's gene order.
Compute the frequency of each k-mer in the chromosome.
Represent each chromosome as a vector in a high-dimensional space where each dimension corresponds to a unique k-mer.
Use dimensionality reduction (PCA, t-SNE, UMAP) to visualize trends.
2. Distance-based Embeddings
Compute pairwise distances between gene orders using a metric like Kendall tau distance, breakpoint distance, or adjacency disruption distance.
Use multidimensional scaling (MDS) or t-SNE/UMAP to embed the chromosomes into a lower-dimensional space.
3. Alignment-based Representation
Use multiple sequence alignment-like methods to align gene orders.
Construct a position-specific frequency matrix (like in sequence logos).
Reduce dimensions using PCA or neural network embeddings.
4. Graph-based Embedding
Construct adjacency graphs from gene orders.
Use graph embeddings (e.g., Node2Vec, GraphSAGE) to obtain vector representations of chromosomes.
5. Next Steps
If you want a simple approach, try k-mer frequency vectors + PCA.
If gene orders have evolutionary significance, try distance-based embeddings + MDS.
If your dataset is large, consider graph-based or learning-based embeddings.
"""

def extract_CDS_gene_orders_only(gene_order_list):
    return [item for item in gene_order_list 
            if not (item[0:5]=='lrRNA' or item[0:5]=='srRNA' or item[0:3]=='Trn')]