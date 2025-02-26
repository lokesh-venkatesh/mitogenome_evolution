from pathlib import Path
from mitofuncs.mito import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sn

gene_locations_folder_path = Path("results/gene_locations")

if not gene_locations_folder_path.exists():
    gene_locations_folder_path.mkdir(parents=True)
    print(f"Folder '{gene_locations_folder_path}' created.")
else: 
    print(f"Folder '{gene_locations_folder_path}' already exists.")

all_species_names = list_all_species_names_from_file_path()

for species_name in all_species_names:
    species_gb_filepath = Path(f"genbank_files/{species_name}/{species_name}_mitochondrion.gb")
    genome_sequence = extract_genome_sequence(genbank_filepath=species_gb_filepath)

    fig, ax = plt.subplots(figsize=(12, 6))            
    gene_df = pd.read_csv(f"genbank_files/{species_name}/{species_name}_modified_cleaned_gene_data.tsv", delimiter="\t").squeeze()
    human_gene_list = list((pd.read_csv(f"genbank_files/Homo_sapiens/Homo_sapiens_modified_cleaned_gene_data.tsv", delimiter="\t").squeeze())['Gene'])
    num_genes = len(human_gene_list)
    cmap = plt.colormaps['rainbow_r'].resampled(num_genes) # Change 'rainbow' to another colormap if desired
    gene_colors = {human_gene_list[i]: cmap(i / num_genes) for i in range(num_genes)}

    for idx, row in gene_df.iterrows():
        # Convert gene positions to the window range
        gene_start = row['Start']
        gene_stop = row['Stop']
        gene_name = row['Gene']
        # Draw a shaded bar for the gene position
        if gene_name in gene_colors:
            colour_shade = gene_colors[gene_name]
        else:
            colour_shade = 'gray'
        plt.axvspan(gene_start, gene_stop, color=colour_shade,
                    alpha=0.6)
        text_x = (gene_start + gene_stop) / 2  # Midpoint
        text_y = ax.get_ylim()[1]/2
        ax.text(text_x, text_y, gene_name, ha='center', va='top', fontsize=9, rotation=90, color='black')

    # Set x-axis limits to match the start and stop of the gene positions
    #x_min = gene_df['Start'].min()
    #x_max = gene_df['Stop'].max()
    #ax.set_xlim(x_min, x_max)
    ax.set_xlim(0,len(genome_sequence))

    # Remove y-axis ticks
    ax.yaxis.set_ticks([])

    plt.xlabel(f"Gene location along {species_name} Genome")
    #plt.ylabel("Corresponding Nucleotide Frequency")
    plt.title(f"Order of all genes in the {species_name} Genome")
    #plt.legend(loc='lower right', frameon=True, fontsize=10)
    plt.savefig(f"results/gene_locations/{species_name}_gene_order.png", dpi=300)
    #plt.show()
    plt.close()
    print(species_name)