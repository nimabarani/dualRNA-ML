#!/bin/bash
#SBATCH --time=02:30:00
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=8GB
#SBATCH --output='../logs/seqkit-%A.out'

for exp_dir in /home/nima/scratch/dual_rna/experiments/*
do
    echo $(basename $exp_dir)
    seqkit stats $exp_dir/*/*/raw_reads/*.fastq
    echo ""
done