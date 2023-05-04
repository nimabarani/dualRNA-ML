#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem-per-cpu=2GB
#SBATCH --output='../logs/sequence_extraction-%A.out'

module load StdEnv/2020 bedtools/2.29.2 bedops/2.4.39

for genome_dir in /home/nima/projects/def-lpenacas/nima/newDual/genomes/pathogens/*
do
    awk -F'\t' '$3 ~ /gene/' $genome_dir/genomic.gtf > $genome_dir/genes/genes.gtf
    
    grep -w -f $genome_dir/genes/gene_ids.txt $genome_dir/genes/genes.gtf > $genome_dir/genes/selected_genes.gtf
    
    gtf2bed < $genome_dir/genes/selected_genes.gtf > $genome_dir/genes/selected_genes.bed
    
    bedtools getfasta -fi $genome_dir/genomic.fna \
                      -bed $genome_dir/genes/selected_genes.bed \
                      -nameOnly \
                      -fo $genome_dir/genes/selected_genes.fna

    awk 'BEGIN{RS=">"}{print $1","$2;}' $genome_dir/genes/selected_genes.fna > $genome_dir/genes/selected_genes.csv
    # echo $(wc -l $genome_dir/genes/selected_genes.csv)
done