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
        if item[0:5]=='lrRNA' or item[0:5]=='srRNA' or item[0:3]=='tRN':
            n_RNA += 1
    n_CDS = n_tot-n_RNA
    return [n_tot, n_CDS, n_RNA]

def extract_CDS_gene_orders_only(gene_order_list):
    return [item for item in gene_order_list 
            if not (item[0:5]=='lrRNA' or item[0:5]=='srRNA' or item[0:3]=='Trn')]

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