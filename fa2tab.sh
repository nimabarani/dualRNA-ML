#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12GB
#SBATCH --output='../logs/fx2tab-%A.out'

CSVS="/home/nima/projects/def-lpenacas/nima/newDual/csvs"
for genome_dir in /home/nima/projects/def-lpenacas/nima/newDual/genomes/hosts/*
do
    if [ -e $genome_dir/genes.fna ]
    then
        refseq="$(basename $genome_dir)"
        awk 'BEGIN{RS=">"}{print $1"\t"$2;}' $genome_dir/genes.fna > $CSVS/$refseq.tab
    fi
done