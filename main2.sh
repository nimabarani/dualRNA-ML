#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5GB
#SBATCH --output='../logs/main-%A.out'

module load java/17.0.2 trimmomatic/0.39 gcc/9.3.0 sra-toolkit/3.0.0 star/2.7.9a bowtie2/2.4.4 bedops/2.4.39 sambamba/0.8.0 subread/2.0.3

fasterQ() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3

    if [ $layout = 'paired' ]; then
        fasterq-dump $sra_dir/$sra_id.sra --split-files -O $sra_dir/

    elif [ $layout = 'single' ]; then
        fasterq-dump $sra_dir/$sra_id.sra -O $sra_dir/
    fi
}

trimming() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3

    if [ $layout = 'paired' ]; then
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 16 \
        $sra_dir/${sra_id}_1.fastq $sra_dir/${sra_id}_2.fastq \
        $sra_dir/${sra_id}_paired_1.fastq $sra_dir/${sra_id}_unpaired_1.fastq \
        $sra_dir/${sra_id}_paired_2.fastq $sra_dir/${sra_id}_unpaired_2.fastq \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -phred33

    elif [ $layout = 'single' ]; then
        java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar SE -threads 16 \
        $sra_dir/${sra_id}.fastq $sra_dir/${sra_id}_trimmed.fastq\
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -phred33
    fi
}

host_alignment() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3
    local genome_dir=$4
    local host_dir=$(dirname $genome_dir)

    if [ $layout = 'paired' ]
    then
        STAR --genomeDir $genome_dir/ \
        --runThreadN 16 \
        --readFilesIn $sra_dir/${sra_id}_paired_1.fastq $sra_dir/${sra_id}_paired_2.fastq \
        --outFileNamePrefix $sra_dir/star/star_ \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile $host_dir/genomic.gtf \
        --quantMode GeneCounts

    elif [ $layout = 'single' ]
    then
        STAR --genomeDir $genome_dir/ \
        --runThreadN 16 \
        --readFilesIn $sra_dir/${sra_id}_trimmed.fastq \
        --outFileNamePrefix $sra_dir/star/star_ \
        --outSAMtype BAM SortedByCoordinate \
        --sjdbGTFfile $host_dir/genomic.gtf \
        --quantMode GeneCounts
    fi
}

pathogen_alignment() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3
    local genome_dir=$4

    if [ $layout = 'paired' ]
    then
        bowtie2 -p 16 --very-sensitive -x $genome_dir/index/pathogen -q \
        -1 $sra_dir/${sra_id}_paired_1.fastq \
        -2 $sra_dir/${sra_id}_paired_2.fastq \
        -S $sra_dir/pathogen_align.sam

    elif [ $layout = 'single' ]
    then
        bowtie2 -p 16 --very-sensitive -x $genome_dir/index/pathogen -q \
        $sra_dir/${sra_id}_trimmed.fastq \
        -S $sra_dir/pathogen_align.sam
    fi
}

gene_quantification() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3
    local genome_dir=$4

    if [ $layout = 'paired' ]
    then
        featureCounts -T 16 -p --countReadPairs -s $strandness \
        -t gene \
        -a $genome_dir/genomic.gtf \
        -o $logs_dir/$working_dir/$experiment_type/${srid}_pathogen.txt \
        $sra_dir/pathogen_sorted_name.bam
    
    elif [ $layout = 'single' ]
    then
        featureCounts -T 16 \
        -t gene \
        -a $genome_dir/genomic.gtf \
        -o $logs_dir/$working_dir/$experiment_type/${srid}_pathogen.txt \
        $sra_dir/pathogen_sorted_name.bam
    fi

}

## ---------------------------------- ##
## -------------  MAIN  ------------- ##
## ---------------------------------- ##

json_data=$(cat setup2.json)

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

            sra_id=$(basename $sra_dir)
            echo "## $sra_id"
            # Print the directory name

            if [[ $experiment_type == 'infection' ]]
            then
                echo "### fasterq-dump"
                fasterQ $layout $sra_dir $sra_id
                echo ""
            
                echo "### trimmomatic"
                trimming $layout $sra_dir $sra_id
                echo ""
            
                echo "### STAR"
                host_alignment $layout $sra_dir $sra_id $host_dir
                sambamba sort -t 16 -n $sra_dir/star/star_Aligned.sortedByCoord.out.bam -o $sra_dir/host_sorted_name.bam
                echo ""
            
                echo "### Bowtie2"
                pathogen_alignment $layout $sra_dir $sra_id $pathogen_dir
                sambamba view -t 16 -S $sra_dir/pathogen_align.sam -f bam -o $sra_dir/pathogen_align.bam
                sambamba sort -t 16 -n $sra_dir/pathogen_align.bam -o $sra_dir/pathogen_sorted_name.bam
                echo ""
            fi
        fi
    done
done
# Loop through each item in the parent directory
