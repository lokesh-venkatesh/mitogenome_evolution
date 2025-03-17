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
        if item[0:5]=='lrRNA' or item[0:5]=='srRNA' or item[0:3]=='Trn':
            n_RNA += 1
    n_CDS = n_tot-n_RNA
    return [n_tot, n_CDS, n_RNA]

def extract_CDS_gene_orders_only(gene_order_list):
    return [item for item in gene_order_list 
            if not (item[0:5]=='lrRNA' or item[0:5]=='srRNA' or item[0:3]=='Trn')]

all_feature_vectors = {species: extract_feature_descriptive_vector(all_species_gene_dataframes[species]['Gene'].tolist()) 
                       for species in all_species_names}

all_CDS_gene_orders = {species: extract_CDS_gene_orders_only(all_species_gene_dataframes[species]['Gene'].tolist()) 
                        for species in all_species_names}

# Extracting the components for plotting
species_names = list(all_feature_vectors.keys())
n_tot_values = [all_feature_vectors[species][0] for species in species_names]
n_CDS_values = [all_feature_vectors[species][1] for species in species_names]
n_RNA_values = [all_feature_vectors[species][2] for species in species_names]

# Creating the subplots
fig, axs = plt.subplots(1, 3, figsize=(24, 8))

# Histogram for total genes
axs[0].hist(n_tot_values, bins=range(min(n_tot_values), max(n_tot_values) + 2, 1), color='b', edgecolor='black')
axs[0].set_title('Total Genes Distribution')
axs[0].set_xlabel('Total Genes')
axs[0].set_ylabel('Frequency')
axs[0].set_xlim(10, 50)
axs[0].set_ylim(0, len(all_species_names))

# Histogram for CDS genes
axs[1].hist(n_CDS_values, bins=range(min(n_CDS_values), max(n_CDS_values) + 2, 1), color='g', edgecolor='black')
axs[1].set_title('CDS Genes Distribution')
axs[1].set_xlabel('CDS Genes')
axs[1].set_ylabel('Frequency')
axs[1].set_xlim(10, 50)
axs[1].set_ylim(0, len(all_species_names))

# Histogram for RNA genes
axs[2].hist(n_RNA_values, bins=range(min(n_RNA_values), max(n_RNA_values) + 2, 1), color='r', edgecolor='black')
axs[2].set_title('RNA Genes Distribution')
axs[2].set_xlabel('RNA Genes')
axs[2].set_ylabel('Frequency')
axs[2].set_xlim(10, 50)
axs[2].set_ylim(0, len(all_species_names))

plt.tight_layout()
plt.savefig(gene_trends_folder_path/'histograms_of_gene_counts.png', dpi=300)