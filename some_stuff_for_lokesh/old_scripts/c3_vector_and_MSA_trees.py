from pathlib import Path
from mitofuncs.mito import *

all_species_names = list_all_species_names_from_file_path()

MSA_and_Vector_trees_folder_path = Path("data/MSA_and_kmer_Vector_trees")
if not MSA_and_Vector_trees_folder_path.exists():
    MSA_and_Vector_trees_folder_path.mkdir(parents=True)

for k in [5,6]:
    dict_of_all_kmer_vectors = {}
    for species_name in all_species_names:
        species_gb_filepath = Path(f"genbank_files/{species_name}/{species_name}_mitochondrion.gb")
        genome_sequence = extract_genome_sequence(genbank_filepath=species_gb_filepath)
        nucl_freq = return_genome_nucl_frequencies_as_dict(genome_sequence)
        species_kmer_vector = generate_full_kmer_vector(k, species_name, genome_sequence, nucl_freq, 
                                  f"{species_name}_{k}mer_vector.TSV", 
                                  folder_path="data/whole_genome_kmer_vectors")
        dict_of_all_kmer_vectors[species_name] = species_kmer_vector
    
    vector_correlation_dataframe = generate_correlation_matrix(dictionary_of_vectors=dict_of_all_kmer_vectors,
                                                               r_type="modified_pcc")

    custom_order = ['Apis_mellifera', 'Meloidogyne_javanica', 'Drosophila_melanogaster', 'Aedes_albopictus', 
                    'Romanomermis_culicivorax', 'Rhipicephalus_sanguineus', 'Anopheles_gambiae', 
                    'Caenorhabditis_elegans', 'Locusta_migratoria', 'Onchocerca_volvulus', 'Florometra_serratissima', 
                    'Ascaris_suum', 'Loligo_bleekeri', 'Albinaria_coerulea', 'Homarus_americanus', 
                    'Katharina_tunicata', 'Limulus_polyphemus', 'Speleonectes_tulumensis', 'Didelphis_virginiana', 
                    'Pisaster_ochraceus', 'Asterias_amurensis', 'Artemia_franciscana', 'Platynereis_dumerilii', 
                    'Cucumaria_miniata', 'Ophiopholis_aculeata', 'Mus_musculus', 'Ornithorhynchus_anatinus', 
                    'Branchiostoma_lanceolatum', 'Branchiostoma_floridae', 'Daphnia_pulex', 'Fasciola_hepatica', 
                    'Urechis_caupo', 'Mytilus_edulis', 'Metridium_senile', 'Lumbricus_terrestris', 'Boa_constrictor', 
                    'Parastichopus_californicus', 'Pelomedusa_subrufa', 'Pateria_pectinifera', 'Sphenodon_punctatus', 
                    'Danio_rerio', 'Dinodon_semicarinatus', 'Lithobates_catesbeianus', 'Cepaea_nemoralis', 
                    'Python_bivittatus', 'Strongylocentrotus_purpuratus', 'Gloydius_blomhoffii', 
                    'Fejervarya_limnocharis', 'Xenopus_tropicalis', 'Alligator_mississippiensis', 
                    'Crocodylus_porosus', 'Homo_sapiens', 'Falco_peregrinus', 'Struthio_camelus', 
                    'Smithornis_sharpei', 'Vidua_chalybeata', 'Gallus_gallus', 'Rhea_americana', 'Bipes_biporus', 
                    'Aythya_americana', 'Balanoglossus_carnosus']
    
    """NOTE THAT THIS CUSTOM ORDER GIVEN ABOVE IS IN TERMS OF INCREASING GC CONTENT"""

    reordered_correlation_dataframe = vector_correlation_dataframe.loc[custom_order, custom_order]
    fig, ax = plt.subplots(figsize=(25,25))
    hm = sn.heatmap(reordered_correlation_dataframe,
                    cmap='coolwarm', cbar_kws={"shrink": 1.0},
                    xticklabels=custom_order, yticklabels=custom_order,
                    square=True, vmin=0, vmax=1,
                    annot=True, fmt=".2g")
    fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
    plt.title(f"Whole genome {k}-mer vectors cross correlation")
    plt.xlabel("species")
    plt.ylabel("species")
    plt.tight_layout()
    plt.savefig(f"data/MSA_and_kmer_Vector_trees/all_{k}mer_species_correlation_heatmap.png", dpi=400)
    plt.close()

    vector_distance_dataframe = generate_distance_matrix(dict_of_all_kmer_vectors, r_type="modified_pcc")
    MEGA_distance_df_string = convert_distance_dataframe_to_MEGA11_file_format(vector_distance_dataframe, 
                                                                    f"All {k}mer Vectors of the Current Dataset")
    write_string_to_file(string=MEGA_distance_df_string, 
                            filepath=f"data/MSA_and_kmer_Vector_trees/all_{k}mer_species.MEG")

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
                         filepath=f"data/FASTA_files/all_species_in_current_dataset.fasta")