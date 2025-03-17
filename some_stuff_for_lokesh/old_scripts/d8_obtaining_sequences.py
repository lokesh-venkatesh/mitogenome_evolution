from pathlib import Path
from mitofuncs.mitoevo import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sn
from difflib import SequenceMatcher
import json

W = 200

output_path = Path("results/levenshtein/all_levenshtein_pairs.tsv")
all_levenshtein_pairs_df = pd.read_csv(output_path, delimiter="\t", index_col=0).squeeze()

all_levenshtein_pairs = []

for key, value in all_levenshtein_pairs_df.items():
    all_levenshtein_pairs.append((key, value))

all_sequences = []

print(len(list_all_species_names_from_file_path()))
print(len(all_levenshtein_pairs))
print(len(all_sequences))

for species_pair in all_levenshtein_pairs:
    species1 = species_pair[0]
    genome1 = extract_genome_sequence(genbank_filepath=f"genbank_files/{species1}/{species1}_mitochondrion.gb")
    species1_gene_data = pd.read_csv(f"genbank_files/{species1}/{species1}_modified_cleaned_gene_data.tsv", 
                        delimiter="\t").squeeze()
    list1 = species1_gene_data['Gene'].tolist()

    species2 = species_pair[1]
    genome2 = extract_genome_sequence(genbank_filepath=f"genbank_files/{species2}/{species2}_mitochondrion.gb")
    species2_gene_data = pd.read_csv(f"genbank_files/{species2}/{species2}_modified_cleaned_gene_data.tsv", 
                        delimiter="\t").squeeze()
    list2 = species2_gene_data['Gene'].tolist()

    all_translocation_events = find_translocated_chunks(list1, list2)

    for translocation in all_translocation_events: #get the start location of the first gene in translocation
        if any(translocation == list1[i:i+len(translocation)] for i in range(len(list1) - len(translocation) + 1)):
            L1 = len(genome1)
            genome1 = genome1*3
            start_location = L1 + species1_gene_data.loc[species1_gene_data['Gene'] == translocation[0], 'Start'].values[0]
            stop_location = L1 + species1_gene_data.loc[species1_gene_data['Gene'] == translocation[-1], 'Stop'].values[0]
            sequence_1 = genome1[start_location-W:start_location]
            sequence_2 = genome1[stop_location:stop_location+W]
            all_sequences.append(sequence_1)
            all_sequences.append(sequence_2)
        elif any(translocation[::-1] == list1[i:i+len(translocation[::-1])] for i in range(len(list1) - len(translocation[::-1]) + 1)):
            L1 = len(genome1)
            genome1 = genome1*3
            start_location = L1 + species1_gene_data.loc[species1_gene_data['Gene'] == translocation[0], 'Stop'].values[0]
            stop_location = L1 + species1_gene_data.loc[species1_gene_data['Gene'] == translocation[-1], 'Start'].values[0]
            sequence_1 = genome1[stop_location:stop_location+W]
            sequence_2 = genome1[start_location-W:start_location]
            all_sequences.append(sequence_1)
            all_sequences.append(sequence_2)
        if any(translocation == list2[i:i+len(translocation)] for i in range(len(list2) - len(translocation) + 1)):
            L2 = len(genome2)
            genome2 = genome2*3
            start_location = L2 + species2_gene_data.loc[species2_gene_data['Gene'] == translocation[0], 'Start'].values[0]
            stop_location = L2 + species2_gene_data.loc[species2_gene_data['Gene'] == translocation[-1], 'Stop'].values[0]
            sequence_1 = genome2[start_location-W:start_location]
            sequence_2 = genome2[stop_location:stop_location+W]
            all_sequences.append(sequence_1)
            all_sequences.append(sequence_2)
        elif any(translocation[::-1] == list2[i:i+len(translocation[::-1])] for i in range(len(list2) - len(translocation[::-1]) + 1)):
            L2 = len(genome2)
            genome2 = genome2*3
            start_location = L2 + species2_gene_data.loc[species2_gene_data['Gene'] == translocation[0], 'Stop'].values[0]
            stop_location = L2 + species2_gene_data.loc[species2_gene_data['Gene'] == translocation[-1], 'Start'].values[0]
            sequence_1 = genome2[stop_location:stop_location+W]
            sequence_2 = genome2[start_location-W:start_location]
            all_sequences.append(sequence_1)
            all_sequences.append(sequence_2)

print(len(list_all_species_names_from_file_path()))
print(len(all_levenshtein_pairs))
print(len(all_sequences))

final_correlation_folder_path = Path("results/final_correlation")
if not final_correlation_folder_path.exists():
    final_correlation_folder_path.mkdir(parents=True)

with open(final_correlation_folder_path / f"all_target_sequences_with_window_size_{W}.json", "w") as f:
    json.dump(all_sequences, f)