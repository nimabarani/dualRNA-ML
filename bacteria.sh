#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24GB
#SBATCH --output='../logs/bacteria-%A.out'

module load bowtie2/2.4.4


bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen
