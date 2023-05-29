#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=7
#SBATCH --output='../logs/bowtie_align-%A.out'

module load bowtie2/2.4.4

align() {
    srid=$1
    genome_dir="../genomes/pathogens/GCF_000022165.1"

    sra_dir="/home/nima/scratch/dual_rna/raw_reads/experiment3/infection/$srid"
    echo $sra_dir
    bowtie2 -p 16 --very-sensitive -x $genome_dir/index/pathogen -q \
    -1 $sra_dir/${srid}_1.fastq \
    -2 $sra_dir/${srid}_2.fastq \
    -S $sra_dir/pathogen_align3.sam
}

align SRR10005563
# align SRR15886835
# align SRR15886836
