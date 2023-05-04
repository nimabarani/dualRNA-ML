#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=32GB
#SBATCH --output='../logs/sambamba-%A.out'

module load sambamba/0.8.0

for sam_dir in /home/nima/scratch/dual_rna/raw_reads/experiment*/{infection,host}/*/star
do
    if [[ ! -f $sam_dir/output_sorted_name.bam ]]
    then
        echo $sam_dir
        sambamba sort -t 18 -n $sam_dir/star_Aligned.sortedByCoord.out.bam -o $sam_dir/output_sorted_name.bam
    fi
done

sam_dir=/home/nima/scratch/dual_rna/raw_reads/experiment10/infection/SRR5581952
sambamba view -t 18 -S $sam_dir/pathogen_align.sam -f bam -o $sam_dir/proband_tmp.bam
sambamba sort -t 18 -n $sam_dir/proband_tmp.bam -o $sam_dir/output_sorted_name.bam

sam_dir=/home/nima/scratch/dual_rna/raw_reads/experiment10/infection/SRR5581954
sambamba view -t 18 -S $sam_dir/pathogen_align.sam -f bam -o $sam_dir/proband_tmp.bam
sambamba sort -t 18 -n $sam_dir/proband_tmp.bam -o $sam_dir/output_sorted_name.bam

# done