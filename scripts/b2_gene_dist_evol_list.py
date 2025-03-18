from mitofuncs.mitoevo import *
from mitofuncs.mitosignal import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

gene_orders_distances_dataframe = pd.read_csv(
    "results/gene_order_distances/all_gene_order_distances.tsv", 
    delimiter="\t", index_col=0)

current_dataset_evol_order = arrange_in_order_based_on_distances(gene_orders_distances_dataframe)
for item in current_dataset_evol_order:
    print(item)


current_dataset_evol_order_series = pd.Series(current_dataset_evol_order)
output_path = Path("data/current_dataset_evol_order.tsv")
current_dataset_evol_order_series.to_csv(output_path, sep='\t', header=False)