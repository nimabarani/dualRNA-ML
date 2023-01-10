#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/samtools-%A.out'

module load samtools/1.16.1

samDir="../raw_data/pathogen3"

samtools sort $samDir/output.sam -o $samDir/output_sorted_position.bam
samtools sort -n $samDir/output.sam -o $samDir/output_sorted_name.bam
