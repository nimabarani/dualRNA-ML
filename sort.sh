#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='./logs/sam-%A.out'


module load StdEnv/2020 samtools/1.16.1


samtools sort -n ./STAR/star_Aligned.sortedByCoord.out.bam -o ./STAR/star.bam
