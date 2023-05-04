#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/bowtie_index-%A.out'

module load bowtie2/2.4.4

genomeDir="/home/nima/projects/def-lpenacas/nima/newDual/genomes/pathogens/GCF_000013425.1"

mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen

