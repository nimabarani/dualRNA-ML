#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4GB
#SBATCH --output='../logs/trim-%A.out'

module load java/17.0.2 trimmomatic/0.39

srid="SRR15886837"
dir="/home/nima/scratch/dual_rna/raw_reads/experiment1/host/$srid"

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 6 \
$dir/${srid}_1.fastq $dir/${srid}_2.fastq \
$dir/${srid}_paired_1.fastq $dir/${srid}_unpaired_1.fastq \
$dir/${srid}_paired_2.fastq $dir/${srid}_unpaired_2.fastq \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -phred33
