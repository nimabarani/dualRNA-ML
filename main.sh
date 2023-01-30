#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=24GB
#SBATCH --output='../logs/main-%A.out'

module load java/17.0.2 trimmomatic/0.39 gcc/9.3.0 sra-toolkit/3.0.0 star/2.7.9a bowtie2/2.4.4 bedops/2.4.39

srid=$1
sra_dir="../raw_reads/$srid"
layout=$2
host_dir=$3
pathogen_dir=$4



# Loop through each item in the parent directory
for item in ../*
do
  # Check if the item is a directory
  if [ -d $item ]
  then
    # Print the directory name
    fasterQ
    trimming
    host_alignment
    pathogen_alignment

    echo $item
  fi
done


function fasterQ{
    if [ $layout = 'p' ]; then
        fasterq-dump $sra_dir/$srid.sra --split-files -O $sra_dir/
    elif [ $layout = 's' ]; then
        fasterq-dump $sra_dir/$srid.sra -O $sra_dir/
    fi
}

function trimming {
    if [ $layout = 'p' ]; then
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 6 \
        $sra_dir/${srid}_1.fastq $sra_dir/${srid}_2.fastq \
        $sra_dir/${srid}_paired_1.fastq $sra_dir/${srid}_unpaired_1.fastq \
        $sra_dir/${srid}_paired_2.fastq $sra_dir/${srid}_unpaired_2.fastq \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    elif [ $layout = 's' ]; then
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 6 \
        $sra_dir/${srid}.fastq $sra_dir/${srid}_trimmed.fastq\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    fi
}

function host_alignment {
    if [ $layout = 'p' ]; then
        STAR --genomeDir $host_dir/index/ \
        --runThreadN 6 \
        --readFilesIn $sra_dir/${srid}_paired_1.fastq $sra_dir/${srid}_paired_2.fastq \
        --outFileNamePrefix $sra_dir/star/star_ \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile $host_dir/genomic.gtf \
        --quantMode GeneCounts
    elif [ $layout = 's' ]; then
        STAR --genomeDir $host_dir/index/ \
        --runThreadN 6 \
        --readFilesIn $sra_dir/${srid}_trimmed.fastq \
        --outFileNamePrefix $sra_dir/star/star_ \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile $host_dir/genomic.gtf \
        --quantMode GeneCounts
    fi
}

function pathogen_alignment {
    if [ $layout = 'p' ]; then
        bowtie2 -p 6 --very-sensitive -x $pathogen_dir/index/pathogen -q \
        -1 $sra_dir/${srid}_paired_1.fastq \
        -2 $sra_dir/${srid}_paired_2.fastq \
        -S $sra_dir/pathogen_align.sam
    elif [ $layout = 's' ]; then
        bowtie2 -p 6 --very-sensitive -x $pathogen_dir/index/pathogen -q \
        $sra_dir/${srid}_trimmed.fastq \
        -S $sra_dir/pathogen_align.sam
    fi
}