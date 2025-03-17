from pathlib import Path
from mitofuncs.mitoevo import *

all_species_names = list_all_species_names_from_file_path()

moving_window_vector_analysis_folder_path = Path("data/mv_wndw_vctr_analysis")
if not moving_window_vector_analysis_folder_path.exists():
    moving_window_vector_analysis_folder_path.mkdir(parents=True)

moving_window_vector_analysis_folder_path = Path("results/mv_wndw_vctr_analysis")
if not moving_window_vector_analysis_folder_path.exists():
    moving_window_vector_analysis_folder_path.mkdir(parents=True)

for species_name in all_species_names:
    for k in [5,6]:
        species_gb_filepath = Path(f"genbank_files/{species_name}/{species_name}_mitochondrion.gb")
        genome_sequence = extract_genome_sequence(genbank_filepath=species_gb_filepath)

        L = len(genome_sequence)
        genome_sequence = genome_sequence*2

        window_corresponding_vectors = {}

        for i in range(0,L,100):
            window = genome_sequence[i:i+1000]
            nucl_freq = return_genome_nucl_frequencies_as_dict(window)
            window_corresponding_vectors[i] = generate_reduced_kmer_vector(generate_full_kmer_vector(k, f"{species_name}_{i}", window, nucl_freq, 
                                                f"{species_name}_{k}mer_vector_window_positioned_at_{i}_in_genome.TSV", 
                                                folder_path="data/mv_wndw_vctr_analysis"))
        species_moving_window_vector_analysis = pd.DataFrame.from_dict(window_corresponding_vectors, orient="index").T
        species_moving_window_vector_analysis = species_moving_window_vector_analysis.map(lambda x: np.log10(x) if x > 0 else np.log10(0.001))
        #fig, ax = plt.subplots(figsize=(25,25))
        hm = sn.heatmap(species_moving_window_vector_analysis,
                        cmap='coolwarm' #cbar_kws={"shrink": 1.0},
                        #square=True,
                        #vmin=0, vmax=1,
                        #annot=True, fmt=".2g"
                        )
        #fig.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)
        plt.title(f"{species_name} {k}-mer vectors moving window analysis")
        plt.xlabel("window position index")
        plt.ylabel("motifs along the vector")
        plt.tight_layout()
        plt.savefig(f"results/mv_wndw_vctr_analysis/{species_name}_{k}mer_vector_moving_window_heatmap.png", dpi=400)
        plt.close()
        print(f"{species_name}, {k}mer-vector moving window analysis done")