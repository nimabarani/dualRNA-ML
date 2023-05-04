#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --output='../logs/slurm-%A.out'
#SBATCH --mem-per-cpu 8G

module load  StdEnv/2020  gcc/9.3.0  sra-toolkit/2.10.8

for sra_dir in /home/nima/scratch/dual_rna/raw_reads/experiment14/infection/ERR560444/*
do
    srid=$(basename $sra_dir)
    echo "--- ${srid} ---"
    fasterq-dump $sra_dir/$srid.sra --split-files -O $sra_dir/
done
