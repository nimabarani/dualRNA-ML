#!/bin/bash
#SBATCH --time=02:30:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/star-%A.out'


module load star/2.7.9a

align() {
    srid=$1
    sra_dir="/home/nima/scratch/dual_rna/raw_reads/exp/infection/${srid}"
    genome_dir="/home/nima/projects/def-lpenacas/nima/newDual/genomes/hosts/GCF_000001405.40"
    
    STAR --genomeDir $genome_dir/index-150 \
    --runThreadN 16 \
    --readFilesIn $sra_dir/${srid}_paired_1.fastq $sra_dir/${srid}_paired_2.fastq \
    --outFileNamePrefix $sra_dir/star2/\
    --outReadsUnmapped Fastx \
    --sjdbGTFfile $genome_dir/genomic.gtf

    # STAR --genomeDir $genome_dir/ \
    # --runThreadN 6 \
    # --readFilesIn $sra_dir/${srid}_trimmed.fastq \
    # --outFileNamePrefix $sra_dir/star/star_ \
    # --outSAMtype BAM SortedByCoordinate \
    # --sjdbGTFfile $host_dir/genomic.gtf \
    # --quantMode GeneCounts
}

align SRR15886834
align SRR15886835
align SRR15886836
