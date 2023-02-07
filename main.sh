#!/bin/bash
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=6
#SBATCH --mem-per-cpu=32GB
#SBATCH --output='../logs/main-%A.out'

module load java/17.0.2 trimmomatic/0.39 gcc/9.3.0 sra-toolkit/3.0.0 star/2.7.9a bowtie2/2.4.4 bedops/2.4.39

fasterQ() {
    local layout=$1
    local sra_dir=$2
    local srid=$3

    if [ $layout = 'paired' ]; then
        fasterq-dump $sra_dir/$srid.sra --split-files -O $sra_dir/

    elif [ $layout = 'single' ]; then
        fasterq-dump $sra_dir/$srid.sra -O $sra_dir/
    fi
}

trimming() {
    local layout=$1
    local sra_dir=$2
    local srid=$3

    if [ $layout = 'paired' ]; then
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 6 \
        $sra_dir/${srid}_1.fastq $sra_dir/${srid}_2.fastq \
        $sra_dir/${srid}_paired_1.fastq $sra_dir/${srid}_unpaired_1.fastq \
        $sra_dir/${srid}_paired_2.fastq $sra_dir/${srid}_unpaired_2.fastq \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -phred64

    elif [ $layout = 'single' ]; then
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 6 \
        $sra_dir/${srid}.fastq $sra_dir/${srid}_trimmed.fastq\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -phred64
    fi
}

host_alignment() {
    local layout=$1
    local sra_dir=$2
    local srid=$3
    local genome_dir=$4
    local host_dir=$(dirname $genome_dir)

    if [ $layout = 'paired' ]
    then
        STAR --genomeDir $genome_dir/ \
        --runThreadN 6 \
        --readFilesIn $sra_dir/${srid}_paired_1.fastq $sra_dir/${srid}_paired_2.fastq \
        --outFileNamePrefix $sra_dir/star/star_ \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile $host_dir/genomic.gtf \
        --quantMode GeneCounts

    elif [ $layout = 'single' ]
    then
        STAR --genomeDir $genome_dir/ \
        --runThreadN 6 \
        --readFilesIn $sra_dir/${srid}_trimmed.fastq \
        --outFileNamePrefix $sra_dir/star/star_ \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile $host_dir/genomic.gtf \
        --quantMode GeneCounts
    fi
}

pathogen_alignment() {
    local layout=$1
    local sra_dir=$2
    local srid=$3
    local genome_dir=$4
    # pathogen_dir=$(dirname $genome_dir)

    if [ $layout = 'paired' ]
    then
        bowtie2 -p 6 --very-sensitive -x $genome_dir/index/pathogen -q \
        -1 $sra_dir/${srid}_paired_1.fastq \
        -2 $sra_dir/${srid}_paired_2.fastq \
        -S $sra_dir/pathogen_align.sam

    elif [ $layout = 'single' ]
    then
        bowtie2 -p 6 --very-sensitive -x $genome_dir/index/pathogen -q \
        $sra_dir/${srid}_trimmed.fastq \
        -S $sra_dir/pathogen_align.sam
    fi
}



json_data=$(cat setup.json)


for data in $(echo $json_data | jq -c '.[]')
do
    # Extract the values of each key
    working_dir=$(echo $data | jq -r .working_dir)
    layout=$(echo $data | jq -r '.layout')
    pathogen_genome_dir=$(echo $data | jq -r '.pathogen_genome_dir')
    pathogen_dir=/home/nima/projects/def-lpenacas/nima/newDual/genomes/pathogens/$pathogen_genome_dir

    host_genome_dir=$(echo $data | jq -r '.host_genome_dir')
    host_dir=/home/nima/projects/def-lpenacas/nima/newDual/genomes/hosts/$host_genome_dir

    # Print the values
    echo "------------------------ $working_dir ------------------------"
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
            echo "---------------- SRA: $srid ----------------"
            # Print the directory name

            echo "-------- BEGINED: fasterq-dump --------"
            fasterQ $layout $sra_dir $srid
            echo "-------- FINISHED: fasterq-dump --------"
            echo ""

            echo "-------- BEGINED: Trimmomatic --------"
            trimming $layout $sra_dir $srid
            echo "-------- FINISHED: Trimmomatic --------"
            echo ""

            if [[ $result == 'host' || $result == 'infection' ]]
            then
                echo "-------- BEGINED: STAR --------"
                host_alignment $layout $sra_dir $srid $host_dir
                echo "-------- FINISHED: STAR --------"
                echo ""

            elif [[ $result == 'pathogen' || $result == 'infection' ]]
            then
                echo "-------- BEGINED: Bowtie2 --------"
                pathogen_alignment $layout $sra_dir $srid $pathogen_dir
                echo "-------- FINISHED: Bowtie2 --------"
                echo ""
            fi
        fi
    done
done
# Loop through each item in the parent directory



