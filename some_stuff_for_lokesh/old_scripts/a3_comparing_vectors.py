from pathlib import Path
from mitofuncs.mito import *

folder_path = Path("results/whole_vector_comparisons")

if not folder_path.exists():
    folder_path.mkdir(parents=True)
    print(f"Folder '{folder_path}' created.")
else: 
    print(f"Folder '{folder_path}' already exists.")

all_animal_species = list_all_species_names_from_file_path()

for k in [5,6]:
    all_species_kmer_vectors = {}
    for species in all_animal_species:
        genome = extract_genome_sequence(f"genbank_files/{species}/{species}_mitochondrion.gb")
        nucl_freq = return_genome_nucl_frequencies_as_dict(genome)
        kmer_vector = generate_full_kmer_vector(k=k, species_name=species, genome_sequence=genome, nucl_freq=nucl_freq, file_name=f"{species}_whole_genome_{k}mer_vector.tsv",folder_path="data/whole_genome_kmer_vectors")
        reduced_kmer_vector = generate_reduced_kmer_vector(kmer_vector)
        all_species_kmer_vectors[species] = reduced_kmer_vector
    correlation_matrix = []
    for species_1, vector_1 in all_species_kmer_vectors.items():
        corrln_row = []
        for species_2, vector_2 in all_species_kmer_vectors.items():
            corrln_row.append(pcc(vector_1, vector_2))
        correlation_matrix.append(corrln_row)
    correlation_dataframe = pd.DataFrame(data=correlation_matrix, index=all_animal_species, columns=all_animal_species)
    custom_order = ['Arabidopsis_thaliana', 
                    'Saccharomyces_cerevisiae', 
                    'Drosophila_melanogaster',
                    'Caenorhabditis_elegans',
                    'Xenopus_tropicalis',
                    'Danio_rerio',
                    'Homo_sapiens', 
                    'Mus_musculus']
    reordered_correlation_dataframe = correlation_dataframe.loc[custom_order, custom_order]
    plt.figure(figsize=(8,8))
    hm = sn.heatmap(reordered_correlation_dataframe,
                    cmap='coolwarm',
                    xticklabels=custom_order, yticklabels=custom_order,
                    square=True, vmin=0, vmax=1,
                    annot=True, fmt=".2g")
    plt.title(f"Whole {k}-mer vectors correlation heatmap")
    plt.xlabel("species")
    plt.ylabel("species")
    plt.tight_layout()
    plt.savefig(f"results/whole_vector_comparisons/whole_{k}mer_vector_heatmap.png", dpi=300)

    distance_dataframe = reordered_correlation_dataframe.map(lambda x: -np.log10(x))
    mega_file_string = convert_dataframe_to_MEGA11_file_format(distance_dataframe, f"all {k}-mer whole vectors")
    
    write_string_to_file(mega_file_string, f"results/whole_vector_comparisons/mega11 {k}-mer full vectors string file.MEG")