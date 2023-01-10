#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2GB
#SBATCH --output='../logs/fasterq-%A.out'

module load gcc/9.3.0  sra-toolkit/3.0.0

fasterq-dump ../raw_data/reads/SRR11307870.sra --split-files -O ../raw_data/reads/fastq_files/SRR11307870
