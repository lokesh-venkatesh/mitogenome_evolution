from mitofuncs.mito import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

species_name = "Homo_sapiens"
species_gb_filepath = Path(f"genbank_files/{species_name}/{species_name}_mitochondrion.gb")
genome_sequence = extract_genome_sequence(genbank_filepath=species_gb_filepath)
L = len(genome_sequence)

genome_sequence += genome_sequence

motif_counts = []

import re
def count_variable_motif(sequence: str, motif: str) -> int:
    regex_pattern = motif.replace("N", "[ATGC]")
    return len(re.findall(regex_pattern, sequence))

for i in range(0,L):
    window = genome_sequence[i:i+100]
    #number_of_motifs = window.count("TTAGGG")
    number_of_motifs = count_variable_motif(window, "GNNNNNNNNNNG")
    motif_counts.append(number_of_motifs)

genome_length = len(motif_counts)
x_positions = list(range(len(motif_counts)))

# Plot the bar chart
plt.figure(figsize=(12, 5))
plt.bar(x_positions, motif_counts, width=1, color='royalblue')

# Labels and title
plt.xlabel("Genome Position")
plt.ylabel("Motif Count")
plt.title("Motif Occurrence Across Genome (Sliding Window)")

# Improve visibility
plt.xlim(0, genome_length)
plt.tight_layout()

# Show plot
plt.savefig("results/TFAM_binding_sites_plot.png",dpi=300)