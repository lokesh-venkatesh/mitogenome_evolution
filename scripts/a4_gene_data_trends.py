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
I AM LEAVING THE FOLLOWING NOTE AS WELL AS THIS CODE FILE AS IT IS TO REVISIT LATER AND DO THE APPROPRIATE PLOTS.

1. PLOT (n(tot), n(CDS), n(RNA)), MAYBE WITH A PCA OR SOMETHING LIKE THAT
2. MAYBE THE NUMBER DOESN'T MATTER AS MUCH, BUT THE PROTEIN CODING GENES MIGHT. SO YOU WILL HAVE TO FIGURE OUT SOME WAY TO REPRESENT THE PROTEIN CODING GENES AS A VECTOR OF SORTS.
THEN PLOT THESE TO SHOW SOME SORT OF LINEAR TREND WITH THE PROTEIN GENES' ORDERS (HERE YOU CAN DEFINITELY TRY PCA OR SOMETHING)
3. TRY K-MEANS CLUSTERING ON SOME OF THESE DATAPOINTS OR WHATEVER. 

TLDR: FIND SOME NICE TRENDS TO VISUALISE!
"""

def extract_feature_descriptive_vector(gene_order_list):
    n_tot = len(gene_order_list)
    n_CDS = 0
    n_RNA = 0
    for item in gene_order_list:
        if item.startswith('lrRNA') or item.startswith('srRNA') or item.startswith('tRN'):
            n_RNA += 1
    n_CDS = n_tot-n_RNA
    return [n_tot, n_CDS, n_RNA]

def extract_CDS_gene_orders_only(gene_order_list):
    return [item for item in gene_order_list 
            if item.startswith('lrRNA') or item.startswith('srRNA') or item.startswith('tRN')]

all_feature_vectors = {species: extract_feature_descriptive_vector(all_species_gene_dataframes[species]['Gene'].tolist()) 
                       for species in all_species_names}

# Extract CDS gene orders for coloring
all_CDS_gene_orders = {species: extract_CDS_gene_orders_only(all_species_gene_dataframes[species]['Gene'].tolist()) 
                        for species in all_species_names}

# Assign colors to unique CDS gene orders
unique_CDS_orders = list(set(tuple(order) for order in all_CDS_gene_orders.values()))
color_map = {order: plt.cm.tab20(i / len(unique_CDS_orders)) for i, order in enumerate(unique_CDS_orders)}

# Create a color list for each species based on their CDS gene orders
colors = [color_map[tuple(all_CDS_gene_orders[species])] for species in all_species_names]

# Convert feature vectors to a DataFrame for easier plotting
feature_df = pd.DataFrame.from_dict(all_feature_vectors, orient='index', columns=['n_tot', 'n_CDS', 'n_RNA'])

# Add jitter to the points to avoid overlap
def add_jitter(arr, scale=0.01):
    return arr + np.random.normal(scale=scale, size=arr.shape)

# Plotting
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# Plot n_tot vs n_CDS
axes[0].scatter(add_jitter(feature_df['n_tot'].values), add_jitter(feature_df['n_CDS'].values), c=colors)
axes[0].set_xlabel('Total number of genes (n_tot)')
axes[0].set_ylabel('Number of protein coding genes (n_CDS)')
axes[0].set_title('n_tot vs n_CDS')

# Plot n_tot vs n_RNA
axes[1].scatter(add_jitter(feature_df['n_tot'].values), add_jitter(feature_df['n_RNA'].values), c=colors)
axes[1].set_xlabel('Total number of genes (n_tot)')
axes[1].set_ylabel('Number of RNA genes (n_RNA)')
axes[1].set_title('n_tot vs n_RNA')

# Plot n_CDS vs n_RNA
axes[2].scatter(add_jitter(feature_df['n_CDS'].values), add_jitter(feature_df['n_RNA'].values), c=colors)
axes[2].set_xlabel('Number of protein coding genes (n_CDS)')
axes[2].set_ylabel('Number of RNA genes (n_RNA)')
axes[2].set_title('n_CDS vs n_RNA')

plt.tight_layout()
plt.savefig("results/gene_trends/trends_in_number_of_CDS_RNA_genes.png", dpi=400)
plt.close()

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

for k in [5,6]:
    dict_of_all_kmer_vectors = {}
    for species_name in all_species_names:
        species_gb_filepath = Path(f"data/genbank_files/{species_name}/{species_name}_mitochondrion.gb")
        genome_sequence = extract_genome_sequence(genbank_filepath=species_gb_filepath)
        nucl_freq = return_genome_nucl_frequencies_as_dict(genome_sequence)
        species_kmer_vector = generate_full_kmer_vector(k, species_name, genome_sequence, nucl_freq, 
                                  f"{species_name}_{k}mer_vector.TSV", 
                                  folder_path="data/dump_for_whole_genome_kmer_vectors")
        dict_of_all_kmer_vectors[species_name] = species_kmer_vector

    # Convert the k-mer vectors dictionary to a DataFrame
    kmer_df = pd.DataFrame.from_dict(dict_of_all_kmer_vectors, orient='index')

    # Perform PCA
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(kmer_df)

    # Create a DataFrame with the PCA results
    pca_df = pd.DataFrame(data=principal_components, columns=['PC1', 'PC2'], index=kmer_df.index)

    # Plot the PCA results
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x='PC1', y='PC2', data=pca_df, hue=pca_df.index, palette='tab20', legend=None)
    plt.title(f'PCA of {k}mer vectors of the current dataset')
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.savefig(f"results/gene_trends/PCA_results_on_{k}mer_vectors.png", dpi=400)
    plt.close()