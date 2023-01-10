#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24GB
#SBATCH --output='../logs/bowtie_index-%A.out'

module load bowtie2/2.4.4

srid="SRR11307870"
genomeDir="../raw_data/pathogen3"
readDir="../raw_data/reads/fastq_files/$srid/${srid}"


bowtie2 -p 4 --very-sensitive -x $genomeDir/index/pathogen -q -1 ${readDir}_paired_1.fastq -2 ${readDir}_paired_2.fastq -S $genomeDir/output.sam
