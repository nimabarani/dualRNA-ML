#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=24GB
#SBATCH --output='../logs/index-%A.out'

genomeDir="../genomes/hosts/human"

module load star/2.7.9a

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir  $genomeDir/75 \
--genomeFastaFiles $genomeDir/genomic.fna \
--sjdbGTFfile $genomeDir/genomic.gtf \
--sjdbOverhang 75
