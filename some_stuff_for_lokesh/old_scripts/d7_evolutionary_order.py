from pathlib import Path
from mitofuncs.mito import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sn
from difflib import SequenceMatcher

"""
WE WILL CONSTRUCT A DICTIONARY OF ALL PAIRS OF GENOMES WITH MINIMSED LEVENSHTEIN DISTANCE, BUT STRICTLY >0
"""

levenshtein_distances_matrix = pd.read_csv("results/levenshtein/all_genomes_levenshtein_distance_df", 
                                           delimiter="\t", index_col=0)
all_species_names = list_all_species_names_from_file_path()
all_levenshtein_pairs = []

for species_1 in all_species_names:
    species_1_distances = levenshtein_distances_matrix.loc[species_1]
    min_non_zero_distance = species_1_distances[species_1_distances > 0].min()
    closest_species = species_1_distances[species_1_distances == min_non_zero_distance].index.tolist()
    for species_2 in closest_species:
        all_levenshtein_pairs.append((species_1, species_2))

unique_levenshtein_pairs = set()
for pair in all_levenshtein_pairs:
    sorted_pair = tuple(sorted(pair))
    unique_levenshtein_pairs.add(sorted_pair)

all_levenshtein_pairs = list(unique_levenshtein_pairs)
output_path = Path("results/levenshtein/all_levenshtein_pairs.tsv")
pd.DataFrame(all_levenshtein_pairs, columns=["Species 1", "Species 2"]).to_csv(output_path, index=False, sep="\t")

species1 = "Homo_sapiens"
species2 = "Smithornis_sharpei"
list1 = pd.read_csv(f"genbank_files/{species1}/{species1}_modified_cleaned_gene_data.tsv", 
                    delimiter="\t").squeeze()['Gene'].tolist()
list2 = pd.read_csv(f"genbank_files/{species2}/{species2}_modified_cleaned_gene_data.tsv", 
                    delimiter="\t").squeeze()['Gene'].tolist()
print(find_translocated_chunks(list1, list2))