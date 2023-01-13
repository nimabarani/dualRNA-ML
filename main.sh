#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=24GB
#SBATCH --output='../logs/main-%A.out'

module load java/17.0.2 trimmomatic/0.39 gcc/9.3.0 sra-toolkit/3.0.0 star/2.7.9a

srid=$1
sra_dir="../raw_reads/$srid"
genome_dir=$3

if [ $2 == 'p' ]; then
    fasterq-dump $sra_dir/$srid.sra --split-files -O $sra_dir/
elif [ $2 == 's' ]; then
    fasterq-dump $sra_dir/$srid.sra -O $sra_dir/
fi
echo 'COMPLETED: FasterQ-dump'

if [ $2 == 'p' ]; then
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 6 \
    $sra_dir/${srid}_1.fastq $sra_dir/${srid}_2.fastq \
    $sra_dir/${srid}_paired_1.fastq $sra_dir/${srid}_unpaired_1.fastq \
    $sra_dir/${srid}_paired_2.fastq $sra_dir/${srid}_unpaired_2.fastq \
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

    echo 'COMPLETED: Trimmomatic'

    STAR --genomeDir $genome_dir/index/ \
    --runThreadN 6 \
    --readFilesIn $sra_dir/${srid}_paired_1.fastq $sra_dir/${srid}_paired_2.fastq \
    --outFileNamePrefix $sra_dir/star/star_ \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbGTFfile $genome_dir/genomic.gtf \
    --quantMode GeneCounts

elif [ $2 == 's' ]; then
    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 6 \
    $sra_dir/${srid}.fastq $sra_dir/${srid}_trimmed.fastq\
    LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    echo 'COMPLETED: Trimmomatic'

    STAR --genomeDir $genome_dir/index/ \
    --runThreadN 6 \
    --readFilesIn $sra_dir/${srid}_trimmed.fastq \
    --outFileNamePrefix $sra_dir/star/star_ \
    --outSAMtype BAM SortedByCoordinate \
    --sjdbGTFfile $genome_dir/genomic.gtf \
    --quantMode GeneCounts
fi


echo 'COMPLETED'


