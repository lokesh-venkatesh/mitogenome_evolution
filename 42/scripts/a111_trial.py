import mitofuncs.mitoevo as mito
from Bio import Entrez
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

print(mito.generate_kmer_motif_list(4))