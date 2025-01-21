#!/bin/bash

#SBATCH --job-name=ECO-ANM            # Job name
#SBATCH -N 1                    # number of nodes
#SBATCH -n 1                    # number of tasks (default: allocates 1 core per task)
#SBATCH -c 16
#SBATCH -t 7-00:00:00           # time in d-hh:mm:s
#SBATCH -p general              # partition
#SBATCH -q public               # QOS
#SBATCH -G a100:1
#SBATCH --mem=40G
#SBATCH -o slurm.%j.out         # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err         # file to save job's STDERR (%j = JobId)
#SBATCH --export=NONE           # Purge the job-submitting shell environment

# Always purge modules to ensure consistent environments
module purge
# Load required modules for job's environment
module load mamba/latest
module load vmd-1.9.4-gcc-12.1.0
module load gcc-12.1.0-gcc-11.2.0
module load mesa-22.0.2-gcc-11.2.0


source activate meld_new

python ANM-highthroughput-pipeline.py

