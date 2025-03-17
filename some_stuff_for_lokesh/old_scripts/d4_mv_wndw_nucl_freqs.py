from pathlib import Path
from mitofuncs.mitoevo import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sn

gene_orders_folder_path = Path("results/gene_orders")

if not gene_orders_folder_path.exists():
    gene_orders_folder_path.mkdir(parents=True)
    print(f"Folder '{gene_orders_folder_path}' created.")
else: 
    print(f"Folder '{gene_orders_folder_path}' already exists.")

all_species_names = list_all_species_names_from_file_path()
freqs = {}

for species_name in all_species_names:
    gene_df = pd.read_csv(f"genbank_files/{species_name}/{species_name}_modified_cleaned_gene_data.tsv", delimiter="\t").squeeze()
    CDS_genes = gene_df[gene_df['Gene'].str.startswith(('ND', 'ATP', 'CO', 'CYTB'))]['Gene'].tolist()
    RNA_genes = gene_df[~gene_df['Gene'].isin(CDS_genes)]['Gene'].tolist()

    human_gene_df = pd.read_csv(f"genbank_files/Homo_sapiens/Homo_sapiens_modified_cleaned_gene_data.tsv", 
                                delimiter="\t").squeeze()
    human_CDS_genes = human_gene_df[human_gene_df['Gene'].str.startswith(('ND', 'ATP', 'CO', 'CYTB'))]['Gene'].tolist()
    human_RNA_genes = human_gene_df[~human_gene_df['Gene'].isin(human_CDS_genes)]['Gene'].tolist()

    num_rna = len(human_RNA_genes)
    cmap_rna = plt.colormaps['rainbow_r'].resampled(num_rna)
    RNA_colors = {human_RNA_genes[i]: cmap_rna(i / num_rna) for i in range(num_rna)}

    num_cds = len(human_CDS_genes)
    cmap_cds = plt.colormaps['rainbow_r'].resampled(num_cds)
    CDS_colors = {human_CDS_genes[i]: cmap_cds(i / num_cds) for i in range(num_cds)}

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw={'hspace': 0.3})

    bar_width_trna = 1 / len(RNA_genes)
    bar_width_cds = 1 / len(CDS_genes)

    for i, gene in enumerate(RNA_genes):
        color = RNA_colors.get(gene, 'grey')
        axes[0].bar(i * bar_width_trna, 1, width=bar_width_trna, color=color, alpha=0.6, align='edge')
        axes[0].text(i * bar_width_trna + bar_width_trna / 2, 0.5, gene, ha='center', va='center', fontsize=8, rotation=90)

    axes[0].set_ylabel("RNA Genes")
    axes[0].set_title(f"RNA Gene Order in {species_name} Genome")
    axes[0].set_xlim(0, 1)
    axes[0].set_ylim(0, 1)
    axes[0].set_xticks([])
    axes[0].set_yticks([])

    for i, gene in enumerate(CDS_genes):
        color = CDS_colors.get(gene, 'grey')
        axes[1].bar(i * bar_width_cds, 1, width=bar_width_cds, color=color, alpha=0.6, align='edge')
        axes[1].text(i * bar_width_cds + bar_width_cds / 2, 0.5, gene, ha='center', va='center', fontsize=8, rotation=90)

    axes[1].set_ylabel("CDS Genes")
    axes[1].set_title(f"CDS Gene Order in {species_name} Genome")
    axes[1].set_xlim(0, 1)
    axes[1].set_xlabel("Genes as they occur from left to right in the Genome")
    axes[1].set_ylim(0, 1)
    axes[1].set_xticks([])
    axes[1].set_yticks([])

    plt.savefig(f"results/gene_orders/{species_name}_RNA_and_CDS_gene_orders.png", dpi=400)
    #plt.show()
    plt.close()