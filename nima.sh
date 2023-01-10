#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --output='../logs/fastx-%A.out'


samDir="../raw_data/host"

fastx_quality_stats -N -i ../raw_data/reads/fastq_files/SRR11307870/SRR11307870_1.fastq -o ./nima.txt
