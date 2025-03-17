from pathlib import Path
from mitofuncs.mitoevo import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sn
import json

W = 200

final_correlation_dump_folder_path = Path("results/final_correlation/final_correlation_dump")
if not final_correlation_dump_folder_path.exists():
    final_correlation_dump_folder_path.mkdir(parents=True)

json_file_path = f'results/final_correlation/all_target_sequences_with_window_size_{W}.json'
with open(json_file_path, 'r') as file:
    target_sequences = json.load(file)

N = len(target_sequences)

for k in [5,6]:
    correlation_matrix = []
    for i in range(N):
        sequence_i = target_sequences[i]
        vector_i = generate_full_kmer_vector(k=k, species_name=f"{i}", 
                                             genome_sequence=sequence_i,
                                             nucl_freq=return_genome_nucl_frequencies_as_dict(sequence_i),
                                             file_name=f"{i}th target sequence",
                                             folder_path=final_correlation_dump_folder_path,
                                             )
        row = []
        for j in range(N):
            sequence_j = target_sequences[j]
            vector_j = generate_full_kmer_vector(k=k, species_name=f"{j}", 
                                             genome_sequence=sequence_j,
                                             nucl_freq=return_genome_nucl_frequencies_as_dict(sequence_j),
                                             file_name=f"{j}th target sequence",
                                             folder_path=final_correlation_dump_folder_path,
                                             )
            r = pearson_correlation_coefficient(vector_i,vector_j)
            row.append(r)
        correlation_matrix.append(row)
    
    correlation_df = pd.DataFrame(correlation_matrix, index=range(N), columns=range(N))
    correlation_df.to_csv(f"results/final_correlation/{k}mer_correlation_df.tsv", sep="\t")
    plt.figure(figsize=(16,16))
    hm = sn.heatmap(correlation_df,
                    cmap='coolwarm',
                    #xticklabels=trial_dataset_names, yticklabels=trial_dataset_names,
                    square=True, 
                    vmin=0, vmax=1,
                    annot=True, fmt=".2g"
                    )
    plt.title(f"Heatmap of all of our target translocation sites' sequences")
    plt.xlabel("Sequence Index")
    plt.ylabel("Sequence Index")
    plt.tight_layout()
    plt.savefig(f"results/final_correlation/{k}mer_correlation_plot_of_all_target_translocation_site_sequences.png", dpi=400)
    plt.close()