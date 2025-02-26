from pathlib import Path
from mitofuncs.mito import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sn

moving_window_nucl_freqs_folder_path = Path("results/moving_window_nucl_freqs")

if not moving_window_nucl_freqs_folder_path.exists():
    moving_window_nucl_freqs_folder_path.mkdir(parents=True)
    print(f"Folder '{moving_window_nucl_freqs_folder_path}' created.")
else: 
    print(f"Folder '{moving_window_nucl_freqs_folder_path}' already exists.")

all_species_names = list_all_species_names_from_file_path()
freqs = {}

for species_name in all_species_names:
    for window_size in [100,1000]:
        window_size_folder = Path(f"results/moving_window_nucl_freqs/window_size_{window_size}")
        if not window_size_folder.exists():
            window_size_folder.mkdir(parents=True)
        species_gb_filepath = Path(f"genbank_files/{species_name}/{species_name}_mitochondrion.gb")
        genome_sequence = extract_genome_sequence(genbank_filepath=species_gb_filepath)

        #window_size=100 #THIS CAME OUT TO BE THE OPTIMUM WINDOW SIZE AFTER PERFORMING THE ELBOW ANALYSIS
        positions, nucleotide_frequencies = sliding_window_nucleotide_frequencies(genome_sequence,
                                                                                window_size=window_size,
                                                                                step_size=1)
        fig, ax = plt.subplots(figsize=(12, 6))            

        gene_df = pd.read_csv(f"genbank_files/{species_name}/{species_name}_modified_cleaned_gene_data.tsv", delimiter="\t").squeeze()

        human_gene_list = list((pd.read_csv(f"genbank_files/Homo_sapiens/Homo_sapiens_modified_cleaned_gene_data.tsv", delimiter="\t").squeeze())['Gene'])
        num_genes = len(human_gene_list)
        cmap = plt.colormaps['rainbow_r'].resampled(num_genes)
         # Change 'rainbow' to another colormap if desired
        gene_colors = {human_gene_list[i]: cmap(i / num_genes) for i in range(num_genes)}
        

        min_y = min(min(nucleotide_frequencies[nuc]) for nuc in nucleotide_frequencies)
        max_y = max(max(nucleotide_frequencies[nuc]) for nuc in nucleotide_frequencies)
        ax.set_ylim(min_y, max_y)
        plt.xlim(0,len(genome_sequence))

        nucleotide_colors = {'A': 'red', 'T': 'royalblue', 'G': 'yellowgreen', 'C': 'orange'}
        for nuc in ['A', 'T', 'G', 'C']:
            plt.plot(positions, nucleotide_frequencies[nuc], 
                    label=f"{nuc}", color=nucleotide_colors[nuc], 
                    linewidth=2, linestyle="--", alpha=0.6)

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
                        alpha=0.4)
            text_x = (gene_start + gene_stop) / 2  # Midpoint
            text_y = max_y - 0.05  # Slightly below top of the plot
            ax.text(text_x, text_y, gene_name, ha='center', va='top', fontsize=9, rotation=90, color='black')

        plt.xlabel(f"Window Position along {species_name} Genome")
        plt.ylabel("Corresponding Nucleotide Frequency")
        plt.title(f"(Smoothened) sliding window analysis of Nucleotide Frequencies in {species_name} Genome (Window={window_size}, Step=1)")
        plt.legend(loc='lower right', frameon=True, fontsize=10)
        plt.grid(True)
        plt.savefig(f"results/moving_window_nucl_freqs/window_size_{window_size}/{species_name}_window_size_{window_size}_line_plots", dpi=300)
        plt.close()
        print(species_name)

'''
species_name = "Homo_sapiens"
species_gb_filepath = Path(f"genbank_files/{species_name}/{species_name}_mitochondrion.gb")
genome_sequence = extract_genome_sequence(genbank_filepath=species_gb_filepath)

window_sizes = [10,20,25,50,100,200,250,500,1000,2000,2500,5000,1000]
variance_values = []

def compute_variance(frequencies):
    """ Compute variance of nucleotide frequencies across genome positions. """
    return {nuc: np.var(frequencies[nuc]) for nuc in frequencies}

for w in window_sizes:
    print(w)
    _, freqs = sliding_window_nucleotide_frequencies(genome_sequence, window_size=w)
    variance = compute_variance(freqs)
    avg_variance = np.mean(list(variance.values()))  # Average variance across A, T, G, C
    variance_values.append(avg_variance)

# Plot variance vs. window size
plt.figure(figsize=(8, 5))
plt.plot(window_sizes, variance_values, marker='o', linestyle='-', color='black')
plt.xlabel("Window Size (bp)")
plt.ylabel("Average Variance in Nucleotide Frequencies")
plt.title("Choosing Optimal Window Size for Nucleotide Analysis")
plt.xscale("log")  # Log scale helps visualize better
plt.grid(True)
plt.savefig(f"results/moving_window_nucl_freqs/window_size_{window_size}/{species_name}_elbow_analysis_for_optimum_window_size", dpi=300)
'''