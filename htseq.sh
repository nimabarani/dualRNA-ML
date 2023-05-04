#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=14GB
#SBATCH --output='../logs/htseq-%A.out'

module load python/3.10.2
source ../env/bin/activate

htseq() {
    sra_id=$1
    genome_dir="/home/nima/projects/def-lpenacas/nima/newDual/scratch/raw_reads/exp/infection/${sra_id}"
    annot_dir="/home/nima/projects/def-lpenacas/nima/newDual/genomes/pathogens/GCF_000010505.1"
    htseq-count -n 4 -f bam -r pos -s yes -c $genome_dir/reads_count.csv \
    $genome_dir/output_sorted_coord.bam \
    $annot_dir/genomic.gtf
}

htseq SRR15886834
htseq SRR15886835
htseq SRR15886836