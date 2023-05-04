#!/bin/bash
#SBATCH --time=00:50:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/trim-%A.out'

# module load fastp/0.23.1
source ../env/bin/activate

trim() {
    sra_id=$1
    dir="/home/nima/scratch/dual_rna/raw_reads/exp/infection/$sra_id"
    # fastp -w 16 -i $dir/${sra_id}_1.fastq -I $dir/${sra_id}_2.fastq -o $sra_dir/${srid}_paired_1.fastq -O $sra_dir/${srid}_paired_2.fastq -h $sra_dir/fastp.html
    cutadapt -q 14 --paired --length 20 --phred33 -o $dir/ $dir/${sra_id}_1.fastq $dir/${sra_id}_2.fastq
    # java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 6 \
    # $dir/${srid}.fastq $dir/${srid}_trimmed.fastq\
    # LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 
}

trim SRR15886834
trim SRR15886835
trim SRR15886836
