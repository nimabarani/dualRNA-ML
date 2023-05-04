#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32GB
#SBATCH --output='../logs/samtools-%A.out'

module load samtools/1.16.1

samDir="/home/nima/scratch/dual_rna/raw_reads/experiment2/host/SRR15827084/star"

samtools sort -@ 6 $samDir/star_Aligned.sortedByCoord.out.bam -o $samDir/output_sorted_position.bam
samtools sort -@ 6 -n $samDir/star_Aligned.sortedByCoord.out.bam -o $samDir/output_sorted_name.bam
