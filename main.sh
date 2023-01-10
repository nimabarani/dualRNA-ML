#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=4GB
#SBATCH --output='../logs/main-%A.out'

module load java/17.0.2 trimmomatic/0.39

srid="SRR11307870"
sra_dir="../raw_reads/$srid"

if [$1 == 'p']
then
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 6 \
    $dir/${srid}_1.fastq $dir/${srid}_2.fastq \
    $dir/${srid}_paired_1.fastq $dir/${srid}_unpaired_1.fastq \
    $dir/${srid}_paired_2.fastq $dir/${srid}_unpaired_2.fastq \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
else
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 6 \
    $dir/${srid}_1.fastq $dir/${srid}_2.fastq \
    $dir/${srid}_paired_1.fastq $dir/${srid}_unpaired_1.fastq \
    $dir/${srid}_paired_2.fastq $dir/${srid}_unpaired_2.fastq \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
fi