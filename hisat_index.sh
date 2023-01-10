#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/hindex-%A.out'

module load StdEnv/2020 hisat2/2.2.1

genomeDir="../raw_data/pathogen"

hisat2-build $genomeDir/pathogen.fna $genomeDir/index/index
