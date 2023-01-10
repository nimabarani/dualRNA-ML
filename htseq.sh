#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8GB
#SBATCH --output='../logs/htseq-%A.out'

module load python/3.10.2
source ../env/bin/activate

genomeDir="../raw_data/host"

htseq-count -f bam -r pos -s reverse \
$genomeDir/output.bam \
$genomeDir/genomic.gtf > $genomeDir/reads_count.txt
