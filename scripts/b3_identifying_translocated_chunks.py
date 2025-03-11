from pathlib import Path
from mitofuncs.mito import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sn

evolutionary_order = pd.read_csv("data/evolutionary_order_series.tsv", delimiter="\t", 
                                 header=None, index_col=0).squeeze().tolist()
gene_data_arranged_in_order = {species: pd.read_csv(f"data/genbank_files/{species}/{species}_cleaned_gene_data.tsv", 
                                              delimiter="\t").squeeze() for species in evolutionary_order}

translocated_chunks_data = []

for i in range(1, len(evolutionary_order)):
    species_1 = evolutionary_order[i-1]
    species_2 = evolutionary_order[i]
    
    gene_order_1 = list(gene_data_arranged_in_order[species_1]['Gene'])
    gene_order_2 = list(gene_data_arranged_in_order[species_2]['Gene'])
    
    translocated_chunks = identify_translocations(gene_order_1, gene_order_2)
    print(translocated_chunks)
    print("")
    translocated_chunks_data.append([i-1, species_1, species_2, translocated_chunks, len(translocated_chunks)])

translocated_chunks_df = pd.DataFrame(translocated_chunks_data, columns=['Index', 'Genome_1', 'Genome_2', 'Translocated_Chunks', 'Count'])
print(translocated_chunks_df)
#translocated_chunks_df.to_csv("data/translocated_chunks.csv", index=False)