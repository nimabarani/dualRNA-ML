#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24GB
#SBATCH --output='../logs/bowtie_index-%A.out'

module load bowtie2/2.4.4

srid="SRR10005567"
genomeDir="../genomes/pathogens/STM"
readDir="/home/nima/scratch/dual_rna/raw_reads/experiment3/infection/$srid"

bowtie2 -p 4 --very-sensitive -x $genomeDir/index/pathogen -q -1 $readDir/${srid}_paired_1.fastq -2 $readDir/${srid}_paired_2.fastq -S $genomeDir/output.sam
# bowtie2 -p 4 --very-sensitive -x $genomeDir/index/pathogen -q $readDir/${srid}_trimmed.fastq -S $readDir/pathogen_align.sam
