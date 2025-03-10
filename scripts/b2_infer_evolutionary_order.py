from pathlib import Path
from mitofuncs.mito import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sn

distance_matrix = pd.read_csv(f"results/DCJ_breakpoint/all_genomes_DCJ_breakpoint_distance_df", delimiter="\t", index_col=0)
genome_names = distance_matrix.columns
distances_numpy_matrix = distance_matrix.to_numpy()
evolutionary_order = infer_evolutionary_order(genome_names, distances_numpy_matrix)

evolutionary_order_series = pd.Series(evolutionary_order, index=range(len(evolutionary_order)))
output_path = Path("data/evolutionary_order_series.tsv")
evolutionary_order_series.to_csv(output_path, header=False, sep="\t")
print(f"Evolutionary order saved to {output_path}")