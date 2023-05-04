#!/bin/bash
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12GB
#SBATCH --output='../logs/test-%A.out'

module load python/3.10.2 bedops/2.4.39

source ../env/bin/activate


json_data=$(cat setup.json)

for data in $(echo $json_data | jq -c '.[]')
do
    # Extract the values of each key
    working_dir=$(echo $data | jq -r .working_dir)
    layout=$(echo $data | jq -r '.layout')

    pathogen_genome_dir=$(echo $data | jq -r '.pathogen_genome_dir')
    pathogen_dir=/home/nima/projects/def-lpenacas/nima/newDual/genomes/pathogens/$pathogen_genome_dir

    # Print the values
    echo "# $working_dir"
    echo "Layout: $layout"
    echo "Host Genome Dir: $pathogen_dir"

    sra_dir="/home/nima/scratch/dual_rna/raw_reads/$working_dir/pathogen"
    first_dir=$(ls -1d $sra_dir/* | head -n 1)

    bam_dir="$first_dir/pathogen_align.sam"

    echo "## $srid"
    echo ""

    if [ -e $bam_dir ]
    then
        infer_experiment.py -r $pathogen_dir/genomic.bed -i $bam_dir
    fi

    echo ""

done
