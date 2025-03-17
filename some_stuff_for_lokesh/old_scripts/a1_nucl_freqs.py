from pathlib import Path
from mitofuncs.mitoevo import *

folder_path = Path("results/nucl_freqs")

if not folder_path.exists():
    folder_path.mkdir(parents=True)
    print(f"Folder '{folder_path}' created.")
else: 
    print(f"Folder '{folder_path}' already exists.")

all_species_names = list_all_species_names_from_file_path()

all_nucleotide_frequencies = {}

for species_name in all_species_names:
    gb_filepath = f"genbank_files/{species_name}/{species_name}_mitochondrion.gb"
    genome_sequence = extract_genome_sequence(genbank_filepath=gb_filepath)
    genome_length = len(genome_sequence)
    nucleotide_frequencies = return_genome_nucl_frequencies_as_dict(genome_sequence=genome_sequence)
    for key in list(nucleotide_frequencies.keys()):
        nucleotide_frequencies[key] = float("%.3f"%(nucleotide_frequencies[key]))
    
    # Colors for the bars
    nucleotide_colors = {'A': 'red', 'T': 'royalblue', 'G': 'yellowgreen', 'C': 'orange'}
    fig, ax = plt.subplots()
    # Create bars with the specified colors
    bars = ax.bar(nucleotide_frequencies.keys(), nucleotide_frequencies.values(), color=[nucleotide_colors[n] for n in nucleotide_frequencies])
    for bar in bars: # Add frequency values on top of the bars
        yval = bar.get_height()
        ax.text(bar.get_x() + bar.get_width()/2, yval + 0.01, f'{yval:.2f}', ha='center', va='bottom')
    ax.set_ylim((0,0.5))
    ax.set_title(f'Nucleotide Frequencies of {species_name} Mt. genome')
    ax.set_xlabel('Nucleotide')
    ax.set_ylabel('Frequency')
    text = f'Genome length: {genome_length}\n GC content: {round(nucleotide_frequencies["G"] + nucleotide_frequencies["C"], 5)}'
    ax.text(3.5, 0.45, text, fontsize=10, ha='right')
    plt.savefig(f"results/nucl_freqs/{species_name}.png", dpi=300)
    plt.close()
    all_nucleotide_frequencies[species_name] = nucleotide_frequencies

nucl_freqs_df = pd.DataFrame.from_dict(all_nucleotide_frequencies, orient="index")
nucl_freqs_df.to_csv("results/nucl_freqs/all_species_freqs.TSV",sep="\t")

# Grouped bar plot
nucleotides = ['A', 'T', 'C', 'G']
x = np.arange(len(nucl_freqs_df))  # the label locations
width = 0.2  # the width of the bars

fig, ax = plt.subplots(figsize=(10, 6))

# Create bars for each nucleotide
for i, nucleotide in enumerate(nucleotides):
    ax.bar(x + i * width, nucl_freqs_df[nucleotide], width, label=nucleotide, color=nucleotide_colors[nucleotide])

# Add labels and title
ax.set_xlabel('Species')
ax.set_ylabel('Frequencies')
ax.set_title('Nucleotide frequencies distributions')
ax.set_xticks(x + width * 1.5)
#ax.set_xticklabels(nucl_freqs_df.index, rotation=45)
ax.set_xticklabels(range(len(nucl_freqs_df.index)))
ax.legend()
plt.ylim((0,0.5))
plt.tight_layout()
plt.savefig(f"results/nucl_freqs/all_freqs_grouped_bar_plot.png", dpi=300)
plt.close()


# Stacked bar plot
fig, ax = plt.subplots(figsize=(10, 6))
# Stack the bars
ax.bar(nucl_freqs_df.index, nucl_freqs_df['A'], label="A", color=nucleotide_colors['A'])
ax.bar(nucl_freqs_df.index, nucl_freqs_df['T'], bottom=nucl_freqs_df['A'], label="T", color=nucleotide_colors['T'])
ax.bar(nucl_freqs_df.index, nucl_freqs_df['C'], bottom=nucl_freqs_df['A'] + nucl_freqs_df['T'], label="C", color=nucleotide_colors['C'])
ax.bar(nucl_freqs_df.index, nucl_freqs_df['G'], bottom=nucl_freqs_df['A'] + nucl_freqs_df['T'] + nucl_freqs_df['C'], label="G", color=nucleotide_colors['G'])
# Add labels and title
ax.set_xlabel('Species')
ax.set_ylabel('Frequencies')
ax.set_title('Nucleotide frequencies distributions')
ax.legend()
ax.set_xticklabels(range(len(nucl_freqs_df.index)))
#plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig(f"results/nucl_freqs/all_freqs_stacked_bar_plot.png", dpi=300)
plt.close()