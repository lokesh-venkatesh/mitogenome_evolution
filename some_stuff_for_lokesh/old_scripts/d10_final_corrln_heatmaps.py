from pathlib import Path
from mitofuncs.mitoevo import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import seaborn as sn
import json

k=5

correlation_df = pd.read_csv(f"results/final_correlation/{k}mer_correlation_df.tsv", delimiter="\t",
                             index_col=0)
correlation_df = correlation_df * 2 - 1
plt.figure(figsize=(16,16))
hm = sn.heatmap(correlation_df,
                cmap='coolwarm',
                #xticklabels=trial_dataset_names, yticklabels=trial_dataset_names,
                square=True, 
                vmin=-1, vmax=1,
                #annot=True, fmt=".2g"
                )
plt.title(f"Heatmap of all of our target translocation sites' sequences")
plt.xlabel("Sequence Index")
plt.ylabel("Sequence Index")
plt.tight_layout()
plt.savefig(f"results/final_correlation/{k}mer_correlation_heatmap_of_all_target_translocation_site_sequences.png", dpi=400)
plt.close()