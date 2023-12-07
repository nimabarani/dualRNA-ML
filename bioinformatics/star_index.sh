#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --output='../logs/index-%A.out'


module load star/2.7.9a

genomeDir="/home/nima/projects/def-lpenacas/nima/newDual/genomes/pathogens/GCF_017498685.1"

STAR --runThreadN 16 \
     --runMode genomeGenerate \
     --genomeDir  $genomeDir/star \
     --genomeFastaFiles $genomeDir/genomic.fna \
     --sjdbGTFfile $genomeDir/genomic.gtf
