from pathlib import Path
from mitofuncs.mito import *

folder_path = Path("results/nucl_freqs")

if not folder_path.exists():
    folder_path.mkdir(parents=True)
    print(f"Folder '{folder_path}' created.")
else: 
    print(f"Folder '{folder_path}' already exists.")

all_species = list_all_species_names_from_file_path()
all_genome_sequences = return_all_mitochondrion_genome_sequences_as_dict(all_species)
all_species_5mer_vectors = {species:generate_full_kmer_vector(species_name=species, k=5, genome_sequence=all_genome_sequences[species], nucl_freq=return_genome_nucl_frequencies_as_dict(all_genome_sequences[species]), file_name=f"{species}_whole_genome_5mer_vector.tsv") for species in all_species}
all_species_6mer_vectors = {species:generate_full_kmer_vector(species_name=species, k=6, genome_sequence=all_genome_sequences[species], nucl_freq=return_genome_nucl_frequencies_as_dict(all_genome_sequences[species]), file_name=f"{species}_whole_genome_6mer_vector.tsv") for species in all_species}