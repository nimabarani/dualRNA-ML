#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24GB
#SBATCH --output='../logs/bowtie_index-%A.out'

module load bowtie2/2.4.4

srid=$1
genomeDir="../genomes/pathogens/sl1344"
readDir="../raw_reads/$srid"


# bowtie2 -p 4 --very-sensitive -x $genomeDir/index/pathogen -q -1 $readDir/${readDir}_paired_1.fastq -2 ${readDir}_paired_2.fastq -S $genomeDir/output.sam
bowtie2 -p 4 --very-sensitive -x $genomeDir/index/pathogen -q $readDir/${srid}_trimmed.fastq -S $readDir/pathogen_align.sam
