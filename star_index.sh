#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=32GB
#SBATCH --output='../logs/index-%A.out'

genomeDir="../genomes/hosts/GCF_000001405.40"

module load star/2.7.9a

STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir  $genomeDir/index-75/ \
--genomeFastaFiles $genomeDir/genomic.fna \
--sjdbGTFfile $genomeDir/genomic.gtf \
--sjdbOverhang 74