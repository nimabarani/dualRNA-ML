#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=1GB
#SBATCH --output='../logs/bed-%A.out'

module load bedops/2.4.39
genome_dir="/home/nima/projects/def-lpenacas/nima/newDual/genomes/hosts/GCF_000001635.27"
if [ -e $genome_dir/selected_genes.gtf ]
then
    gtf2bed < $genome_dir/selected_genes.gtf > $genome_dir/genes.bed
fi
# samtools index your.bam
