#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=40GB
#SBATCH --output='../logs/index-%A.out'

genomeDir="../raw_data/host"

module load StdEnv/2020 star/2.7.9a

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir  $genomeDir/index \
--genomeFastaFiles $genomeDir/genomic.fna \
--sjdbGTFfile $genomeDir/genomic.gtf \
--sjdbOverhang 84
