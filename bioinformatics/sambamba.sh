#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1GB
#SBATCH --output='../logs/sambamba-%A.out'

module load sambamba/0.8.0

sambam() {
    local sra_dir=$1
    local sra_id=$(basename $sra_dir)

    # sambamba view -t 16 -S $sra_dir/mymap.sam -f bam -o $sra_dir/mymap.bam
    sambamba sort -t 16 -n $sra_dir/star2/star_Aligned.out.bam \
                  -o $sra_dir/star2/star_Aligned_name.bam
    # sambamba index -t 18 $sam_dir/output_sorted_coord.bam $sam_dir/pathogen.bam.bai
}


for sra_dir in /home/nima/scratch/dual_rna/raw_reads/experiment2/pathogen/*
do
    sambam $sra_dir
done