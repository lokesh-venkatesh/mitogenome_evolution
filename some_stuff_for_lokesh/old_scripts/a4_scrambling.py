
from pathlib import Path
from mitofuncs.mito import *

folder_path = Path("data/scrambling")
if not folder_path.exists():
    folder_path.mkdir(parents=True)
    print(f"Folder '{folder_path}' created.")
else: 
    print(f"Folder '{folder_path}' already exists.")

for subfolder_name in ["scrambled_genomes", "scramble_vectors", "averaged_vectors_together", "all_avg_vectors"]:
    folder_path = Path(f"data/scrambling/{subfolder_name}")
    if not folder_path.exists():
        folder_path.mkdir(parents=True)
        print(f"Folder '{folder_path}' created.")
    else: 
        print(f"Folder '{folder_path}' already exists.")

all_animal_species = list_all_species_names_from_file_path()
all_genomes = {}
for species in all_animal_species:
    all_genomes[species] = extract_genome_sequence(f"genbank_files/{species}/{species}_mitochondrion.gb")
    nucl_freq = return_genome_nucl_frequencies_as_dict(all_genomes[species])

N = 20 #number of scrambling trials

"""NOTE REDO AND CONTINUE FROM EVERYTHING DOWN BELOW"""

def ith_time_randomisation(species, genome_seq, i, N, folder_path="data/scrambling/scrambled_genomes"):
    """RANDOMISES THE GENOME FOR THE iTH TIME AND SAVES IT AS OUTPUT"""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_path = os.path.join(folder_path, f"{species}_{i}th_scramble.txt") # Define the path to the output file
    if os.path.exists(file_path): # Check if the output file exists
        output = open(file_path).read().strip()
        print(f"Loaded scrambled genome for {i}th iteration from {file_path}.")
    else: # If the file does not exist, generate and save the output
        output = randomise(genome_seq)
        write_string_to_file(output, file_path)
        print(f"Generated and saved output for {i}th scrambled genome.")
    return output

def ith_species_scramble_vector_generation(species, scrambled_seq, nucl_freqs, k, i, N, folder_path="data/scrambling/scramble_vectors"):
    """TAKES THE iTH SCRAMBLED GENOME OF THAT SPECIES AND CONSTRUCTS THE K-MER VECTOR AND SAVES IT"""
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    file_path = os.path.join(folder_path, f"{species}_{i}th_scramble_full_{k}mer_vector.TSV") # Define the path to the output file
    if os.path.exists(file_path): # Check if the output file exists
        output = pd.read_csv(file_path, delimiter="\t", index_col=0)
        print(f"Loaded scrambled genome for {i}th iteration from {file_path}.")
    else: # If the file does not exist, generate and save the output
        output = generate_full_kmer_vector(k, species, scrambled_seq, nucl_freqs, folder_path="data/scrambling/scramble_vectors", file_name=f"{species}_{i}th_scramble_full_{k}mer_vector.TSV")
        print(f"Generated and saved output for {i}th scrambled {k}-mer vector.")
    return output

for k in [5,6]:
    all_avg_kmer_vectors_list = []
    for species in all_animal_species:
        raw_genome = all_genomes[species]
        nucl_frequency_dist = return_genome_nucl_frequencies_as_dict(raw_genome)
        dummy_list = []
        for i in range(1,N+1):
            scrambled_genome = ith_time_randomisation(species, raw_genome, i, N)
            dummy = ith_species_scramble_vector_generation(species, scrambled_genome, nucl_frequency_dist, k, i, N)
            dummy_list.append(dummy)
        avg_species_kmer_vector = sum(dummy_list)/len(dummy_list)
        avg_species_kmer_vector.name = f"{species}"
        avg_species_kmer_vector.to_csv(f"data/scrambling/all_avg_vectors/avg_{species}_{k}mer_vector.TSV",sep="\t")
        all_avg_kmer_vectors_list.append(avg_species_kmer_vector)
    all_avg_kmer_vectors_dataframe = pd.concat(all_avg_kmer_vectors_list, axis=1)
    all_avg_kmer_vectors_dataframe.to_csv(f"data/scrambling/averaged_vectors_together/all_avg_{k}mer_vectors.TSV",sep="\t")