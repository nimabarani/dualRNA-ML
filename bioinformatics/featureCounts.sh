#!/bin/bash
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=1GB
#SBATCH --output='../logs/count-%A.out'

module load subread/2.0.3

data=$1

logs_dir="/home/nima/scratch/dual_rna/count_tables"

# Extract the values of each key
working_dir=$(echo $data | jq -r .working_dir)
layout=$(echo $data | jq -r '.layout')
strandness=$(echo $data | jq -r '.strandness')

pathogen_genome_dir=$(echo $data | jq -r '.pathogen_genome_dir')
pathogen_dir=/home/nima/projects/def-lpenacas/nima/newDual/genomes/pathogens/$pathogen_genome_dir

host_genome_dir=$(echo $data | jq -r '.host_genome_dir')
host_dir=/home/nima/projects/def-lpenacas/nima/newDual/genomes/hosts/$host_genome_dir

# Print the values
echo "# $working_dir"
echo ""
echo "Layout: $layout"
echo "Pathogen Genome Dir: $pathogen_genome_dir"
echo "Host Genome Dir: $host_genome_dir"
echo ""

for sra_dir in ~/scratch/dual_rna/experiments/$working_dir/*/*
do

    # Check if the item is a directory
    if [ -d $sra_dir ]
    then

        experiment_type="${sra_dir%/*}"
        experiment_type="${experiment_type##*/}"

        sra_id=$(basename $sra_dir)
        echo "## $sra_id"
        echo "$logs_dir/$working_dir/$experiment_type/"
        if [ -f $sra_dir/sam-bam/host_sorted_name.bam ]
        then
            # Remove index size (like index-75) from host genome directory
            annot_dir="${host_dir%/*}"
            if [ $layout = 'paired' ]
            then
                featureCounts -T 16 -p --countReadPairs -s $strandness \
                -t gene \
                -a $annot_dir/genomic.gtf \
                -o $logs_dir/$working_dir/$experiment_type/${sra_id}_host.txt \
                $sra_dir/sam-bam/host_sorted_name.bam

            elif [ $layout = 'single' ]
            then
                featureCounts -T 16 \
                -t gene \
                -a $annot_dir/genomic.gtf \
                -o $logs_dir/$working_dir/$experiment_type/${sra_id}_host.txt \
                $sra_dir/sam-bam/host_sorted_name.bam
            fi
        fi

        if [ -f $sra_dir/sam-bam/pathogen_sorted_name.bam ]
        then
            annot_dir=$pathogen_dir
            if [ $layout = 'paired' ]
            then
                featureCounts -T 16 -p --countReadPairs -s $strandness \
                -t gene \
                -a $annot_dir/genomic.gtf \
                -o $logs_dir/$working_dir/$experiment_type/${sra_id}_pathogen.txt \
                $sra_dir/sam-bam/pathogen_sorted_name.bam

            elif [ $layout = 'single' ]
            then
                featureCounts -T 16 \
                -t gene \
                -a $annot_dir/genomic.gtf \
                -o $logs_dir/$working_dir/$experiment_type/${sra_id}_pathogen.txt \
                $sra_dir/sam-bam/pathogen_sorted_name.bam
            fi
        fi
    fi
done
