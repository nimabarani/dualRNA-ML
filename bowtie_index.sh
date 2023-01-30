#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/bowtie_index-%A.out'

module load bowtie2/2.4.4

genomeDir="../genomes/pathogens/sl1344"

bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen
