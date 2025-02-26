from pathlib import Path
from mitofuncs.mito import *

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn

whole_genome_nucl_freqs_folder_path = Path("results/whole_genome_nucl_freqs")

if not whole_genome_nucl_freqs_folder_path.exists():
    whole_genome_nucl_freqs_folder_path.mkdir(parents=True)
    print(f"Folder '{whole_genome_nucl_freqs_folder_path}' created.")
else: 
    print(f"Folder '{whole_genome_nucl_freqs_folder_path}' already exists.")

all_species_names = list_all_species_names_from_file_path()
freqs = {}

i=0
for species_name in all_species_names:
    i+=1
    print(i, species_name)
    species_gb_filepath = Path(f"genbank_files/{species_name}/{species_name}_mitochondrion.gb")
    genome_sequence = extract_genome_sequence(genbank_filepath=species_gb_filepath)
    genome_length = len(genome_sequence)
    nucleotide_frequencies = return_genome_nucl_frequencies_as_dict(genome_sequence=genome_sequence)
    for key in list(nucleotide_frequencies.keys()):
        nucleotide_frequencies[key] = float("%.3f"%(nucleotide_frequencies[key]))
    
    nucleotide_colors = {'A': 'red', 'T': 'royalblue', 'G': 'yellowgreen', 'C': 'orange'}
    fig, ax = plt.subplots()
    bars = ax.bar(nucleotide_frequencies.keys(), nucleotide_frequencies.values(), color=[nucleotide_colors[n] for n in nucleotide_frequencies])
    for bar in bars:
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval + 0.01, f'{yval:.2f}', ha='center', va='bottom')
    ax.set_ylim((0,0.5))
    ax.set_title(f'Nucleotide Frequencies of {species_name} Mt. genome')
    ax.set_xlabel('Nucleotide')
    ax.set_ylabel('Frequency')
    text = f'Genome length: {genome_length}\n GC content: {round(nucleotide_frequencies["G"] + nucleotide_frequencies["C"], 5)}'
    ax.text(3.5, 0.45, text, fontsize=10, ha='right')
    plt.savefig(f"results/whole_genome_nucl_freqs/{species_name}.png", dpi=300)
    plt.close()
    freqs[species_name] = nucleotide_frequencies

whole_genome_nucl_freqs_df = pd.DataFrame.from_dict(freqs, orient="index")
whole_genome_nucl_freqs_df.to_csv("results/whole_genome_nucl_freqs/all_species_freqs.TSV",sep="\t")

freqs = dict(sorted(freqs.items(), key=lambda x: x[1]['G'] + x[1]['C']))
print(list(freqs.keys()))

#---------------------------------------------------------------------------------------

species = list(freqs.keys())
nucleotides = ['A', 'T', 'G', 'C']
data = {n: [freqs[sp][n] for sp in species] for n in nucleotides}
ind = np.arange(len(species))
width = 0.2 
plt.figure(figsize=(15, 5))
for i, nucleotide in enumerate(nucleotides):
    plt.bar(ind + i*width, data[nucleotide], width, label=nucleotide)
plt.xticks(ind + width*1.5, species, rotation='vertical')
plt.ylabel("Nucleotide Frequency")
plt.title("Nucleotide Frequency Distributions of Animalia, arranged in increasing GC content")
plt.legend()
plt.tight_layout()
plt.savefig(f"results/whole_genome_nucl_freqs/all_nucl_freqs_bar_plot_distrbns.png", dpi=400)

#---------------------------------------------------------------------------------------

plt.figure(figsize=(15, 5))
bottom = np.zeros(len(species))
for nucleotide in nucleotides:
    values = [freqs[sp][nucleotide] for sp in species]
    plt.bar(species, values, bottom=bottom, label=nucleotide)
    bottom += values
plt.xticks(rotation=90)
plt.ylabel("Nucleotide Frequency")
plt.title("Nucleotide Frequency Stacks of Animalia, arranged in increasing GC content")
plt.legend()
plt.tight_layout()
plt.savefig(f"results/whole_genome_nucl_freqs/all_nucl_freqs_stacks.png", dpi=400)

#---------------------------------------------------------------------------------------
'''
def create_radar_chart(values, title):
    categories = list(values.keys())
    N = len(categories)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    values_list = list(values.values())
    values_list += values_list[:1]
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(polar=True))
    ax.plot(angles, values_list, 'o-', linewidth=2)
    ax.fill(angles, values_list, alpha=0.25)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels(categories)
    ax.set_title(title)
    plt.show()

create_radar_chart(freqs['Homo_sapiens'], 'Homo_sapiens Nucleotide Frequencies (Radar Chart)')

#---------------------------------------------------------------------------------------

data_list = []
for sp, freqs_dict in freqs.items():
    for nucleotide, freq in freqs_dict.items():
        data_list.append({'species': sp, 'nucleotide': nucleotide, 'frequency': freq})
df = pd.DataFrame(data_list)
plt.figure(figsize=(8, 6))
sn.boxplot(x='nucleotide', y='frequency', data=df)
plt.title("Distribution of Nucleotide Frequencies Across Species")
plt.ylabel("Frequency")
plt.show()

#---------------------------------------------------------------------------------------

df_wide = pd.DataFrame.from_dict(freqs, orient='index')
df_wide.index.name = 'species'
df_wide.reset_index(inplace=True)
sn.pairplot(df_wide, vars=['A', 'T', 'G', 'C'])
plt.suptitle("Pairplot of Nucleotide Frequencies", y=1.02)
plt.show()
'''
#---------------------------------------------------------------------------------------

df_heat = pd.DataFrame.from_dict(freqs, orient='index')
plt.figure(figsize=(10, 12))
sn.heatmap(df_heat, annot=True, cmap="coolwarm_r")
plt.title("Heatmap of Nucleotide Frequency Distributions across species")
plt.xlabel("Nucleotide")
plt.ylabel("Species")
plt.tight_layout()
plt.savefig(f"results/whole_genome_nucl_freqs/all_nucl_freqs_heatmap.png", dpi=400)