#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/ebi-%A.out'

# for file in /home/nima/scratch/dual_rna/experiment20/*.gz
# do
#     gzip -d $file
# done

base_dir=/home/nima/scratch/dual_rna/raw_reads/experiment20

cat $base_dir/ERR3419048.fastq $base_dir/ERR3419051.fastq > SRR3419048.fastq
cat $base_dir/ERR3419049.fastq $base_dir/ERR3419052.fastq > SRR3419049.fastq
cat $base_dir/ERR3419050.fastq $base_dir/ERR3419053.fastq > SRR3419050.fastq
