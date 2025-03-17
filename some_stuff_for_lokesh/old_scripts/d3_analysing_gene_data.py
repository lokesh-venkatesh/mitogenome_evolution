from mitofuncs.mitoevo import *
import pandas as pd
from pathlib import Path
from Bio import SeqIO
import re
#from Bio import pairwise2
#from Bio.pairwise2 import format_alignment
from Bio.Align import PairwiseAligner

all_species = list_all_species_names_from_file_path()

'''
product_to_gene_mapping = pd.read_csv("results/cleaned_gene_data/0_product_to_gene_mapping.tsv", sep='\t', header=0).squeeze()
mapping_series = product_to_gene_mapping.set_index("Product")["Gene"].to_dict()
gene_names = list(mapping_series.values())
all_species_gene_data_details = pd.DataFrame(columns=gene_names)
all_species_gene_data_details = all_species_gene_data_details.loc[:, ~all_species_gene_data_details.columns.duplicated()]

for species in all_species:
    species_gene_data = pd.read_csv(Path(f"genbank_files/{species}/{species}_cleaned_gene_data.tsv"), delimiter="\t")
    gene_count = {gene: species_gene_data[species_gene_data['Gene'] == gene].shape[0] for gene in gene_names}
    gene_count['species'] = species
    all_species_gene_data_details = pd.concat([all_species_gene_data_details, pd.DataFrame([gene_count]).set_index('species')], ignore_index=False)
all_species_gene_data_details['total'] = all_species_gene_data_details.sum(axis=1)
all_species_gene_data_details.to_csv(Path("data/all_species_gene_data_details.tsv"), sep="\t")

average_gene_counts = all_species_gene_data_details.mean(axis=0).tolist()
average_gene_counts = [round(count) for count in average_gene_counts]
mismatch_rows = all_species_gene_data_details[all_species_gene_data_details.apply(lambda row: not all(row == average_gene_counts), axis=1)]
mismatch_rows.to_csv(Path("data/anomalous_gene_details.tsv"), sep="\t")

def extract_trna_anticodons(genbank_filepath):
    """Similar to my genbank parse function but exclusively for extracting anticodon sequences"""
    record = SeqIO.read(genbank_filepath, "genbank")
    genome_seq = str(record.seq)  
    genome_length = len(genome_seq)
    data = []
    for feature in record.features:
        if feature.type in ["tRNA"]:
            product = feature.qualifiers.get("product", [None])[0]
            gene = feature.qualifiers.get("gene", [None])[0]
            start = int(feature.location.start) + 1
            stop = int(feature.location.end)
            strand = "+" if feature.location.strand == 1 else "-"
            if "codon_recognized" in feature.qualifiers:
                codon_recognized = feature.qualifiers.get("codon_recognized", [None])[0]
                if codon_recognized.endswith(')'):
                    codon_recognized = codon_recognized[0:3].upper()
                codon_recognized = ''.join(['N' if nucleotide not in 'AUTGC' else nucleotide for nucleotide in codon_recognized])
                data.append([gene, product, start, stop, strand, codon_recognized])
            elif "anticodon" in feature.qualifiers:
                anticodon = feature.qualifiers["anticodon"][0].split("seq:")[-1].strip('"')
                codon_recognized = anticodon.replace('A', 'U').replace('T', 'A').replace('G', 'C').replace('C', 'G')[::1]
                if codon_recognized.endswith(')'):
                    codon_recognized = codon_recognized[0:3].upper()
                codon_recognized = ''.join(['N' if nucleotide not in 'AUTGC' else nucleotide for nucleotide in codon_recognized])
                codon_recognized = codon_recognized.replace('A', 'U').replace('T', 'A').replace('G', 'C').replace('C', 'G')[::-1]
                data.append([gene, product, start, stop, strand, codon_recognized])
            elif gene is not None:
                if gene[3:5] == "L1":
                    codon_recognized = "UUN"
                    data.append([gene, product, start, stop, strand, codon_recognized])
                elif gene[3:5] == "L2":
                    codon_recognized = "CUN"
                    data.append([gene, product, start, stop, strand, codon_recognized])
                elif gene[3:5] == "S1":
                    codon_recognized = "UCN"
                    data.append([gene, product, start, stop, strand, codon_recognized])
                elif gene[3:5] == "S2":
                    codon_recognized = "AGN"
                    data.append([gene, product, start, stop, strand, codon_recognized])
            elif gene is not None:
                if gene[-1] == ')':
                    anticodon = gene[-4:-1].upper()
                    anticodon = anticodon.replace('A', 'U').replace('T', 'A').replace('G', 'C').replace('C', 'G')[::-1]
                    if anticodon.endswith(')'):
                        anticodon = anticodon[0:3].upper()
                    codon_recognized = ''.join(['N' if nucleotide not in 'AUTGC' else nucleotide for nucleotide in anticodon])
                    #codon_recognized = codon_recognized.replace('A', 'U').replace('T', 'A').replace('G', 'C').replace('C', 'G')[::-1]
                    data.append([gene, product, start, stop, strand, codon_recognized])
            else:
                codon_recognized = None
                data.append([gene, product, start, stop, strand, codon_recognized])
    df = pd.DataFrame(data, columns=["Gene", "Product", "Start", "Stop", "Strand", "Codon_recognized"
                                     ])
    df = df[df['Product'].isin(['tRNA-Leu', 'tRNA-Ser'])]
    df = df.drop_duplicates()
    #df = df.fillna(' ')
    df = df.sort_values(by='Codon_recognized', na_position='last').drop_duplicates(subset=["Gene", "Product", "Start", "Stop", "Strand"], keep='first')
    return df


if Path("data/all_species_trna_details.txt").exists():
    Path("data/all_species_trna_details.txt").unlink()
for species in all_species:
    genbank_filepath = Path(f"genbank_files/{species}/{species}_mitochondrion.gb")
    trna_deets = extract_trna_anticodons(genbank_filepath)
    with open("data/all_species_trna_details.txt", "a") as f:
        f.write(f"Species: {species}\n")
        f.write(trna_deets.to_string(index=False))
        f.write("\n\n")
    trna_deets.to_csv(Path(f"genbank_files/{species}/{species}_tRNA_codon_refs.tsv"), sep="\t", index=False)
    species_gene_data = pd.read_csv(Path(f"genbank_files/{species}/{species}_cleaned_gene_data.tsv"), delimiter="\t")
    for _, trna_row in trna_deets.iterrows():
        if pd.notnull(trna_row['Codon_recognized']):
            matching_rows = species_gene_data[species_gene_data['Start'] == trna_row['Start']]
            for idx in matching_rows.index:
                species_gene_data.at[idx, 'Gene'] = f"{species_gene_data.at[idx, 'Gene']} ({trna_row['Codon_recognized']})"
    species_gene_data.to_csv(Path(f"genbank_files/{species}/{species}_modified_cleaned_gene_data.tsv"), sep="\t", index=False)



all_species = list_all_species_names_from_file_path()

product_to_gene_mapping = pd.read_csv("results/cleaned_gene_data/0_product_to_gene_mapping.tsv", sep='\t', header=0).squeeze()
mapping_series = product_to_gene_mapping.set_index("Product")["Gene"].to_dict()
gene_names = list(mapping_series.values())
all_species_gene_data_details = pd.DataFrame(columns=gene_names)
all_species_gene_data_details = all_species_gene_data_details.loc[:, ~all_species_gene_data_details.columns.duplicated()]

for species in all_species:
    species_gene_data = pd.read_csv(Path(f"genbank_files/{species}/{species}_cleaned_gene_data.tsv"), delimiter="\t")
    gene_count = {gene: species_gene_data[species_gene_data['Gene'] == gene].shape[0] for gene in gene_names}
    gene_count['species'] = species
    all_species_gene_data_details = pd.concat([all_species_gene_data_details, pd.DataFrame([gene_count]).set_index('species')], ignore_index=False)
all_species_gene_data_details['total'] = all_species_gene_data_details.sum(axis=1)
all_species_gene_data_details.to_csv(Path("data/all_species_gene_data_details.tsv"), sep="\t")

average_gene_counts = all_species_gene_data_details.mean(axis=0).tolist()
average_gene_counts = [round(count) for count in average_gene_counts]
mismatch_rows = all_species_gene_data_details[all_species_gene_data_details.apply(lambda row: not all(row == average_gene_counts), axis=1)]
mismatch_rows.to_csv(Path("data/anomalous_gene_details.tsv"), sep="\t")


TrnS_UCN_Sequences = {}
TrnS_AGN_Sequences = {}
TrnL_CUN_Sequences = {}
TrnL_UUN_Sequences = {}

for species in all_species:
    species_gene_data = pd.read_csv(Path(f"genbank_files/{species}/{species}_modified_cleaned_gene_data.tsv"), delimiter="\t")
    species_genome = extract_genome_sequence(genbank_filepath=Path(f"genbank_files/{species}/{species}_mitochondrion.gb"))
    if 'TrnS (UCN)' in species_gene_data['Gene'].values:
        rows = species_gene_data[species_gene_data['Gene'] == 'TrnS (UCN)']
        starts = rows['Start'].values
        stops = rows['Stop'].values
        for start, stop in zip(starts, stops):
            subsequence = species_genome[start-1:stop]
            TrnS_UCN_Sequences[species] = (subsequence)
    if 'TrnS (AGN)' in species_gene_data['Gene'].values:
        rows = species_gene_data[species_gene_data['Gene'] == 'TrnS (AGN)']
        starts = rows['Start'].values
        stops = rows['Stop'].values
        for start, stop in zip(starts, stops):
            subsequence = species_genome[start-1:stop]
            TrnS_AGN_Sequences[species] = (subsequence)
    if 'TrnL (CUN)' in species_gene_data['Gene'].values:
        rows = species_gene_data[species_gene_data['Gene'] == 'TrnL (CUN)']
        starts = rows['Start'].values
        stops = rows['Stop'].values
        for start, stop in zip(starts, stops):
            subsequence = species_genome[start-1:stop]
            TrnL_CUN_Sequences[species] = subsequence
    if 'TrnL (UUN)' in species_gene_data['Gene'].values:
        rows = species_gene_data[species_gene_data['Gene'] == 'TrnL (UUN)']
        starts = rows['Start'].values
        stops = rows['Stop'].values
        for start, stop in zip(starts, stops):
            subsequence = species_genome[start-1:stop]
            TrnL_UUN_Sequences[species] = subsequence


def align_sequences(seq1, seq2):
    aligner = PairwiseAligner()
    alignments = aligner.align(seq1, seq2)
    best_alignment = alignments[0]
    score = best_alignment.score
    normalized_score = score / max(len(seq1), len(seq2))
    return normalized_score

trnS_datalist = []
for seq in TrnS_AGN_Sequences.values():
    trnS_datalist.append(seq)
for seq in TrnS_UCN_Sequences.values():
    trnS_datalist.append(seq)
trnS_dataset = {i:trnS_datalist[i] for i in range(len(trnS_datalist))}

trnL_datalist = []
for seq in TrnL_CUN_Sequences.values():
    trnL_datalist.append(seq)
for seq in TrnL_UUN_Sequences.values():
    trnL_datalist.append(seq)
trnL_dataset = {i:trnL_datalist[i] for i in range(len(trnL_datalist))}


query_species = "Cepaea_nemoralis"
query_species_gene_df = pd.read_csv(Path(f"genbank_files/{query_species}/{query_species}_modified_cleaned_gene_data.tsv"), delimiter="\t")
query_species_genome = extract_genome_sequence(Path(f"genbank_files/{query_species}/{query_species}_mitochondrion.gb"))

query_sequences = []
gene = 'TrnL'

gene_rows = query_species_gene_df[query_species_gene_df['Gene'].str.startswith(gene)]
for _, row in gene_rows.iterrows():
    start = row['Start']
    stop = row['Stop']
    sequence = query_species_genome[start-1:stop]
    query_sequences.append([row['Gene'], start, stop, sequence])

og_dataset = trnS_dataset if gene=='TrnS' else trnL_dataset
trial_dataset = {}
for key, value in og_dataset.items():
    trial_dataset[key] = value
for i in range(len(query_sequences)):
    trial_dataset[f"{query_sequences[i][0]}_{str(query_sequences[i][1])}"] = query_sequences[i][3]

trial_dataset_names = list(trial_dataset.keys())
trial_dataset_seqs = list(trial_dataset.values())
correlation_matrix = []
for i in range(len(trial_dataset_seqs)):
    row = []
    for j in range(len(trial_dataset_seqs)):
        row.append(align_sequences(trial_dataset_seqs[i], trial_dataset_seqs[j]))
    correlation_matrix.append(row)

correlation_df = pd.DataFrame(correlation_matrix, index=trial_dataset_names, columns=trial_dataset_seqs)
plt.figure(figsize=(8,8))
hm = sn.heatmap(correlation_df,
                cmap='coolwarm',
                xticklabels=trial_dataset_names, yticklabels=trial_dataset_names,
                square=True, 
                #vmin=0, vmax=1,
                #annot=True, fmt=".2g"
                )
plt.title(f"Clustermap of {query_species} unknown {gene} sequences")
plt.xlabel("seq")
plt.ylabel("seq")
plt.tight_layout()
plt.savefig(f"data/tRNA_editing/{query_species}_{gene}_sequences.png", dpi=400)
#plt.show()
plt.close()

for seq in query_sequences:
    print(seq)
'''

for species in all_species:
    gene_data_filepath = f"genbank_files/{species}/{species}_modified_cleaned_gene_data.tsv"
    gene_data_df = pd.read_csv(gene_data_filepath, delimiter="\t")
    with open("data/all_species_finally_modified_gene_details.txt", "a") as f:
        f.write(f"Species: {species}\n")
        f.write(gene_data_df.to_string(index=False))
        f.write("\n\n")