#!/bin/bash
#SBATCH --time=02:30:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24GB
#SBATCH --output='../logs/bowtie_align-%A.out'

module load bowtie2/2.4.4

align() {
    srid=$1
    genome_dir="../genomes/pathogens/GCF_000010505.1"

    sra_dir="/home/nima/scratch/dual_rna/raw_reads/exp/infection/$srid"
    echo $sra_dir
    bowtie2 -p 16 --very-sensitive -x $genome_dir/index/pathogen -q \
    -1 $sra_dir/${srid}_1.fastq \
    -2 $sra_dir/${srid}_2.fastq \
    -S $sra_dir/pathogen_align3.sam
}

align SRR15886834
# align SRR15886835
# align SRR15886836
