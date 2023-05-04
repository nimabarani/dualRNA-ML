#!/bin/bash
#SBATCH --time=00:10:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/index-%A.out'


module load star/2.7.9a

genomeDir="/home/nima/projects/def-lpenacas/nima/newDual/genomes/pathogens/GCF_000010505.1"

STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir  $genomeDir/index2 \
--genomeFastaFiles $genomeDir/genomic.fna \
--genomeSAindexNbases 9 \
--sjdbGTFfile $genomeDir/genomic.gtf