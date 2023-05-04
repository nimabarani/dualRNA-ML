#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=4
#SBATCH --output='../logs/count-%A.out'

module load subread/2.0.3

json_data=$(cat setup2.json)
logs_dir="/home/nima/scratch/dual_rna/counts"

for data in $(echo $json_data | jq -c '.[]')
do
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

    for sra_dir in ~/scratch/dual_rna/raw_reads/$working_dir/*/*
    do

        # Check if the item is a directory
        if [ -d $sra_dir ]
        then

            experiment_type="${sra_dir%/*}"
            experiment_type="${experiment_type##*/}"

            srid=$(basename $sra_dir)
            echo "## $srid"
            echo "$logs_dir/$working_dir/$experiment_type/"
            if [ -f $sra_dir/host_sorted_name.bam ]
            then
                annot_dir="${host_dir%/*}"
                
                if [ $layout = 'paired' ]
                then
                    featureCounts -T 16 -p --countReadPairs -s $strandness \
                    -t gene \
                    -a $annot_dir/genomic.gtf \
                    -o $logs_dir/$working_dir/$experiment_type/${srid}_host.txt \
                    $sra_dir/host_sorted_name.bam

                elif [ $layout = 'single' ]
                then
                    featureCounts -T 16 \
                    -t gene \
                    -a $annot_dir/genomic.gtf \
                    -o $logs_dir/$working_dir/$experiment_type/${srid}_host.txt \
                    $sra_dir/host_sorted_name.bam
                fi
            fi

            if [[ -f $sra_dir/pathogen_sorted_name.bam ]]
            then
                if [ $layout = 'paired' ]
                then
                    featureCounts -T 16 -p --countReadPairs -s $strandness \
                    -t gene \
                    -a $pathogen_dir/genomic.gtf \
                    -o $logs_dir/$working_dir/$experiment_type/${srid}_pathogen.txt \
                    $sra_dir/pathogen_sorted_name.bam
                
                elif [ $layout = 'single' ]
                then
                    featureCounts -T 16 \
                    -t gene \
                    -a $pathogen_dir/genomic.gtf \
                    -o $logs_dir/$working_dir/$experiment_type/${srid}_pathogen.txt \
                    $sra_dir/pathogen_sorted_name.bam
                fi
            fi
        fi
    done
done
