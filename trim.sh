#!/bin/bash
#SBATCH --time=01:50:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1GB
#SBATCH --output='../logs/trim-%A.out'

# module load fastp/0.23.1
# source ../env/bin/activate
module load java/17.0.2 trimmomatic/0.39

trimming() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3

    if [ $layout = 'paired' ]; then
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 16 \
        $sra_dir/${sra_id}_1.fastq $sra_dir/${sra_id}_2.fastq \
        $sra_dir/${sra_id}_paired_1.fastq $sra_dir/${sra_id}_unpaired_1.fastq \
        $sra_dir/${sra_id}_paired_2.fastq $sra_dir/${sra_id}_unpaired_2.fastq \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -phred33

    elif [ $layout = 'single' ]; then
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 16 \
        $sra_dir/${sra_id}.fastq $sra_dir/${sra_id}_trimmed.fastq\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -phred33
    fi
}
for sra_dir in /home/nima/scratch/dual_rna/raw_reads/exp/{infection,pathogen}/*
do
    sra_id=$(basename $sra_dir)
    trimming paired $sra_dir $sra_id
done
