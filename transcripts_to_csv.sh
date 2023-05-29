#!/bin/bash

# Activate the virtual environment
source /home/nima/projects/def-lpenacas/nima/newDual/env/bin/activate

for input_file in /home/nima/projects/def-lpenacas/nima/newDual/genomes/pathogens/*
do
    echo $(basename $input_file)
    # Remove parenthesis around the gene names and only keep the gene names
    # perl -pe 's/>.*\((.*)\).*/>$1/' $input_file/rna.fna > $input_file/genes/modified_rna.fna

    # Make gene sequence in one line since NCBI RNA files contains
    # 80 characters in each line
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),substr($0, 2));N++;next;} {printf("%s",$0);} END {printf("\n");}' $input_file/genes/modified_rna.fna > $input_file/genes/modified_rna.csv

    echo "python script"
    python gene_selection.py -i $input_file/genes/modified_rna.csv \
                             -g $input_file/genes/gene_ids.csv \
                             -o $input_file/genes
    
    echo "tab2fx"
    # seqkit tab2fx $input_file/genes/selected_genes.csv | seqkit seq -u > $input_file/genes/selected_genes.fna
done