#!/bin/bash
#SBATCH --partition=gpu             # Request GPU partition
#SBATCH --job-name=all_genomes_dwnld         # Name of the job
#SBATCH --output=jobs/%J.out        # Standard output log file (Job ID specific)
#SBATCH --error=jobs/%J.err         # Standard error log file (Job ID specific)
#SBATCH --time=8-00:00:00           # Maximum time (8 days)
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task=8           # Number of CPU cores per task
#SBATCH --mem=4G                    # Memory per node

# Load necessary modules (if required)
# module load python/3.9  # Uncomment if your HPC requires this

# Activate the virtual environment
source /home/madhu/miniconda3/envs/mito_env

# Define paths
SCRIPT_PATH="$1"  # First argument to the script is the Python script path
OUTPUT_DIR="/storage/madhu/lokesh/42"

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Extract script filename (without path)
SCRIPT_NAME=$(basename "$SCRIPT_PATH")

# Run Python script and pass the output directory dynamically
python3 "$SCRIPT_PATH" --output_dir "$OUTPUT_DIR"

# Deactivate virtual environment
conda deactivate