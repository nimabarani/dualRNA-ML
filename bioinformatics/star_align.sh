#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=8GB
#SBATCH --output='../logs/star-%A.out'


module load star/2.7.9a

alignment() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3
    local genome_dir=$4
    local host_dir=$(dirname $genome_dir)

    if [ $layout = 'paired' ]
    then
        STAR --genomeDir $genome_dir/ \
        --runThreadN 16 \
        --readFilesIn $sra_dir/${sra_id}_paired_1.fastq $sra_dir/${sra_id}_paired_2.fastq \
        --outFileNamePrefix $sra_dir/star4/star_ \
        --outSAMtype BAM Unsorted \
        --sjdbGTFfile $host_dir/genomic.gtf

    elif [ $layout = 'single' ]
    then
        STAR --genomeDir $genome_dir/ \
        --runThreadN 16 \
        --readFilesIn $sra_dir/${sra_id}_trimmed.fastq \
        --outFileNamePrefix $sra_dir/star4/star_ \
        --outSAMtype BAM Unsorted \
        --sjdbGTFfile $host_dir/genomic.gtf
    fi
}

genome_dir="/home/nima/projects/def-lpenacas/nima/newDual/genomes/hosts/GCF_001660625.2/index-150"

for sra_dir in /home/nima/scratch/dual_rna/raw_reads/experiment2/infection/*
do
    sra_id=$(basename $sra_dir)
    alignment paired $sra_dir $sra_id $genome_dir
done