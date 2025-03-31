#!/bin/bash
#SBATCH --partition=gpu             # Request GPU partition
#SBATCH --job-name=lv_tst_jb          # Name of the job
#SBATCH --output=/storage/madhu/lokesh/42_outputs/%J.out  # Standard output log file (Job ID specific)
#SBATCH --error=/storage/madhu/lokesh/42_outputs/%J.err   # Standard error log file (Job ID specific)
#SBATCH --time=8-00:00:00           # Maximum time (8 days)
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task=8           # Number of CPU cores per task
#SBATCH --mem=1G                    # Memory per node

echo "test_job_1 run" > /storage/madhu/lokesh/42_outputs/output_test.txt    # Print a message and system information
hostname >> /storage/madhu/lokesh/42_outputs/output_test.txt
date >> /storage/madhu/lokesh/42_outputs/output_test.txt