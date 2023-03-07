#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=32GB
#SBATCH --output='../logs/sambamba-%A.out'

module load sambamba/0.8.0


samDir="/home/nima/scratch/dual_rna/raw_reads/experiment2/host/SRR15827083/star"

sambamba sort -t 12 -n $samDir/star_Aligned.sortedByCoord.out.bam -o $samDir/output_sorted_name.bam 