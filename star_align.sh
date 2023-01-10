#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=40GB
#SBATCH --output='../logs/star-%A.out'


module load star/2.7.9a

srid="SRR11307870"
fastqDir="../raw_data/reads/fastq_files/$srid/${srid}"
hostDir="../raw_data/host"

STAR --genomeDir $hostDir/index/ \
--runThreadN 6 \
--readFilesIn ${fastqDir}_paired_1.fastq ${fastqDir}_paired_2.fastq \
--outFileNamePrefix $hostDir/star2/star_ \
--outSAMtype BAM SortedByCoordinate \
--sjdbGTFfile $hostDir/genomic.gtf \
--quantMode GeneCounts

