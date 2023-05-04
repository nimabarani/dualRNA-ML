#!/bin/bash

for input_file in /home/nima/projects/def-lpenacas/nima/newDual/genomes/hosts/*
do

    # Remove parenthesis around the gene names and only keep the gene names
    perl -pe 's/>.*\((.*)\).*/>$1/' $input_file/rna.fna > $input_file/modified_rna.fna

    # Make gene sequence in one line since NCBI RNA files contains
    # 80 characters in each line
    awk '/^>/ {printf("%s%s\t",(N>0?"\n":""),$0);N++;next;} {printf("%s",$0);} END {printf("\n");}' $input_file/modified_rna.fna > $input_file/output.csv

    # Remove > from the begining of FASTA file
    perl -pe 's/>//' $input_file/output.csv > $input_file/genes.csv

    # Convert tab to comma as a delimiter
    perl -pe 's/\t/,/' $input_file/genes.csv > $input_file/human_genes.csv

done