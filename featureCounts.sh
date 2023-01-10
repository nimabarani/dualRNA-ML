#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/count-%A.out'

module load subread/2.0.3

genomeDir="../raw_data/pathogen3"

featureCounts -T 4 -p --countReadPairs -s 2 \
  -a $genomeDir/genomic.gtf \
  -o $genomeDir/featurecounts.txt \
  $genomeDir/output.sorted_name

