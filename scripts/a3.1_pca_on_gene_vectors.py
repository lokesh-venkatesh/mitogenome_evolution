"""
THIS FILE IS MEANT TO GENERATE PLOTS (VERY ROUGH IDEAS) FOR THE 
VARIOUS POSISBLE TRENDS IN THE GENE COUNTS AND GENE ORDERS
"""

from mitofuncs.mito import *
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
5. Learning-Based Embedding
Train a neural network (e.g., autoencoder) to learn a latent representation of gene orders.
This requires defining a meaningful loss function based on gene order similarity.
6. Next Steps
If you want a simple approach, try k-mer frequency vectors + PCA.
If gene orders have evolutionary significance, try distance-based embeddings + MDS.
If your dataset is large, consider graph-based or learning-based embeddings.
"""

def extract_CDS_gene_orders_only(gene_order_list):
    return [item for item in gene_order_list 
            if not (item[0:5]=='lrRNA' or item[0:5]=='srRNA' or item[0:3]=='Trn')]

for k in [5, 6]:
    dict_of_all_kmer_vectors = {}
    gene_order_labels = []

    for species_name in all_species_names:
        species_gb_filepath = Path(f"data/genbank_files/{species_name}/{species_name}_mitochondrion.gb")
        genome_sequence = extract_genome_sequence(genbank_filepath=species_gb_filepath)
        nucl_freq = return_genome_nucl_frequencies_as_dict(genome_sequence)
        species_kmer_vector = generate_full_kmer_vector(k, species_name, genome_sequence, nucl_freq, 
                                  f"{species_name}_{k}mer_vector.TSV", 
                                  folder_path="data/dump_for_whole_genome_kmer_vectors")
        dict_of_all_kmer_vectors[species_name] = species_kmer_vector

        # Extract CDS gene orders and use them as labels
        gene_order_list = extract_CDS_gene_orders_only(all_species_gene_dataframes[species_name]['Gene'].to_list())
        print(gene_order_list)
        gene_order_labels.append('-'.join(gene_order_list))

    # Convert the k-mer vectors dictionary to a DataFrame
    kmer_df = pd.DataFrame.from_dict(dict_of_all_kmer_vectors, orient='index')

    # Perform PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(kmer_df)

    # Create a DataFrame with the PCA results
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=kmer_df.index)
    pca_df['GeneOrder'] = gene_order_labels

    # Plot the PCA results
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='PC1', y='PC2', data=pca_df, hue='GeneOrder', palette='tab20', legend=None)
    plt.title(f'PCA of {k}mer vectors of the current dataset')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.savefig(f"results/gene_trends/PCA_results_on_{k}mer_vectors_colored_by_gene_order.png", dpi=400)
    plt.close()