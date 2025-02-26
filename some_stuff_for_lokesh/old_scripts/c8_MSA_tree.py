from pathlib import Path
from mitofuncs.mito import *

all_species_names = list_all_species_names_from_file_path()

all_species_FASTA_file = ""

for species_name in all_species_names:
    species_gb_filepath = Path(f"genbank_files/{species_name}/{species_name}_mitochondrion.gb")
    genome_sequence = extract_genome_sequence(genbank_filepath=species_gb_filepath)
    FASTA_files_folder_path = Path("data/FASTA_files")
    if not FASTA_files_folder_path.exists():
        FASTA_files_folder_path.mkdir(parents=True)
    species_FASTA_file_string = f">{species_name}\n{genome_sequence}"
    all_species_FASTA_file += species_FASTA_file_string+"\n"
    write_string_to_file(string=species_FASTA_file_string,
                         filepath=f"data/FASTA_files/{species_name}.fasta")

write_string_to_file(string=all_species_FASTA_file,
                         filepath=f"summary/trees/all_genomes_in_current_dataset.FASTA")