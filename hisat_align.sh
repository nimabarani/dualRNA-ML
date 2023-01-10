#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/hisat-align-%A.out'

module load StdEnv/2020 hisat2/2.2.1

genomeDir="../raw_data/pathogen"
srid="SRR6493960"
readDir="../raw_data/reads/fastq_files/$srid"

hisat2 -x $genomeDir/index/index  $readDir/${srid}_Unmapped.fastq -S $genomeDir/alignment.sam
