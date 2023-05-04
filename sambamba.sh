#!/bin/bash
#SBATCH --time=9:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24GB
#SBATCH --output='../logs/sambamba-%A.out'

module load sambamba/0.8.0

# sambam() {
#     sra_id=$1
#     sam_dir="/home/nima/projects/def-lpenacas/nima/newDual/scratch/raw_reads/experiment2/infection/${sra_id}"
#     sambamba view -t 16 -S $sam_dir/pathogen_align.sam -f bam -o $sam_dir/pathogen_align.bam
#     sambamba sort -t 16 -n $sam_dir/pathogen_align.bam -o $sam_dir/pathogen_sorted_name.bam
#     # sambamba index -t 18 $sam_dir/output_sorted_coord.bam $sam_dir/pathogen.bam.bai
# }
# sambam SRR15886834
# sambam SRR15886835
# sambam SRR15886836


for dir in /home/nima/projects/def-lpenacas/nima/newDual/scratch/raw_reads/experiment{2,4}/*/*/pathogen_align.sam
do
    sam_dir=$(dirname $dir)
    sambamba view -t 16 -S $sam_dir/pathogen_align.sam -f bam -o $sam_dir/pathogen_align.bam
    sambamba sort -t 16 -n $sam_dir/pathogen_align.bam -o $sam_dir/pathogen_sorted_name.bam
done