#!/bin/bash
#SBATCH --time=10:30:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=32GB
#SBATCH --output='../logs/sambamba-%A.out'

module load sambamba/0.8.0

for sam_dir in /home/nima/scratch/dual_rna/raw_reads/experiment*/{infection,pathogen}/*/*.sam
do
    sra_dir=$(dirname $sam_dir)

    sambamba view -t 18 -S $sra_dir/pathogen_align.sam -f bam -o $sra_dir/proband_tmp.bam
    sambamba sort -t 18 -n $sra_dir/proband_tmp.bam -o $sra_dir/output_sorted_name.bam

done 
# sambamba view -t 12 -S $sam_dir/pathogen_align.sam -f bam -o $sam_dir/proband_tmp.bam
# sambamba sort -t 12 -n $sam_dir/proband_tmp.bam -o $sam_dir/output_sorted_name.bam
# done