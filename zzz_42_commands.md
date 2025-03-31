## Connecting to 42:
`ssh your_username@192.168.1.152` 

## Navigating directories:
To open a directory: `cd /storage/madhu/` 
To print the current directory: `pwd`
Make a new directory: `mkdir folder_name`
List the contents of your current directory:
`ls`
`ls -1`
Open a specific file in a local editor: `cat output_test.txt`

## Using Slurm
To list current running jobs using *slurm*: `squeue`
To list current jobs under your username: `squeue -u madhu`
To check for a specific job ID: `squeue -j 99535`

## Creating python scripts
To create a new .py file: `nano my_script.py`
(you can use anyone of nano, vim or vi)
To save the file, do Ctrl+X, then Y, then Enter (this is for Nano)

## Creating a bash script
To create a new .sh file: `nano job.slurm`

---
#### Contents to put inside the bash script:
#!/bin/bash
#SBATCH --job-name=my_python_job
#SBATCH --output=output/job_output.log
#SBATCH --error=output/job_error.log
#SBATCH --time=00:10:00  # Adjust as needed (HH:MM:SS)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G  # Adjust as needed

module load python/3.8  # Load Python module (if necessary)... Adjust based on available modules
mkdir -p output

python my_script.py
---

## Submitting jobs:
To submit the job: `sbatch test_job.slurm` (it will return a job ID - this can also be seen on squeue)
Checking job status: `squeue -u madhu` 
To check a specific job's details: `scontrol show job 12345` (where 12345 is the corresponding job ID)

## Checking job output:
Navigating to the output directory once the job is complete: `ls output/`
(You should see three things, job_output.log, job_error.log, result.txt)

To view the contents of a .txt file: `cat output/result.txt`
(If job fails, then check the error log, `cat test_job.err`)

## Transfering files and results between local desktop and HPC:
Copying results from 42 to desktop: 
`scp -r madhu@192.168.1.152:/storage/madhu/lokesh/42_outputs C:\Users\lokes\Desktop\mitoevo`

HPC to Laptop: `scp your_username@192.168.1.152:/storage/madhu/lk_test.slurm ~/Desktop/`
Laptop to HPC: `scp ~/Desktop/my_script.sh your_username@192.168.1.152:/storage/madhu/lokesh/`
Copying an entire directory: `scp -r C:\Users\lokes\Desktop\mitoevo madhu@192.168.1.152:/storage/madhu/lokesh/`

## Other stuff about the 42 cluster:
Listing available nodes:
`cat /etc/hosts`
`cat /etc/cluster/config`
If it's a slurm cluster, list available nodes using: `sinfo`

## Check disk space (overall storage):
`df -h`

## Removing files from directory:
*Permanently delete* a specific file within a subfolder: `rm /storage/madhu/lokesh/script.py`
*Permanently delete multiple files* at once: `rm /storage/madhu/lokesh/file1.py /storage/madhu/lokesh/file2.py`
(listing each of them one by one)
**DANGEROUS**... *permanently delete entire directory*: `rm -r /storage/madhu/lokesh/old_project` (-r recursively deletes everything in the directory)
Use -i in front to confirm each deletion, basically interactive mode: `rm -ri /storage/madhu/lokesh/old_project`

Force remove a protected file: `rm -f /home/madhu/lokesh/script.py`
Force remove a directory: `rm -rf /home/madhu/lokesh/old_project`

---
## Updated bash script for usage that Aparajita gave me:
#!/bin/bash
#SBATCH --partition=gpu             # Request GPU partition
#SBATCH --job-name=test_job         # Name of the job
#SBATCH --output=jobs/%J.out        # Standard output log file (Job ID specific)
#SBATCH --error=jobs/%J.err         # Standard error log file (Job ID specific)
#SBATCH --time=8-00:00:00           # Maximum time (8 days)
#SBATCH --ntasks=1                  # Number of tasks
#SBATCH --cpus-per-task=8           # Number of CPU cores per task
#SBATCH --mem=1G                    # Memory per node

echo "Hello from Slurm!" > output_test.txt    #Print a message and system information
hostname >> output_test.txt
date >> output_test.txt
---

## Checking list of modules on a HPC:
`module avail`
`pip list`
To check libraries using conda: `conda list`
To check list of environments: `conda env list`

# NOTE: all of the first nine lines in the bash script remain the same when you're submitting stuff. It's the stuff below that will change depending on what you want. For example: 
`# Cancel the job with Job ID 99535`
`scancel 99535`
`# Confirm deletion (optional)`
`echo "Job 99535 has been deleted." > delete_job_log.txt`