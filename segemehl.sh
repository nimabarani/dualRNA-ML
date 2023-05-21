#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=5GB
#SBATCH --output='../logs/segemehl-%A.out'

module load sambamba/0.8.0

align() {
    local sra_dir=$1
    local genome_dir=$2
    sra_id=$(basename $sra_dir)

    echo $sra_dir
    segemehl.x -t 16 \
               -i $genome_dir/genomic.idx \
               -d $genome_dir/genomic.fna \
               -q $sra_dir/raw_reads/ERR5530739.fastq \
               -o $sra_dir/mymap.sam &> $sra_dir/logs/segemehl.out

    # sambamba view -t 16 -S $sra_dir/mymap.sam -f bam -o $sra_dir/mymap.bam
    # sambamba sort -t 16 -n $sra_dir/mymap.bam -o $sra_dir/mymap_sorted.bam
}
genome_dir="../genomes/pathogens/GCF_000013425.1"

# for sra_dir in /home/nima/scratch/dual_rna/raw_reads/experiment2/infection/SRR15827087
# do
#     align $sra_dir $genome_dir
# done
sra_dir="/home/nima/scratch/dual_rna/experiments/experiment9/pathogen/ERR5530739"
align $sra_dir $genome_dir