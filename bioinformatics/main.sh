#!/bin/bash
#SBATCH --time=00:40:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=24GB
#SBATCH --output='../logs/main-%A.out'

module load java/17.0.2 gcc/9.3.0 star/2.7.9a bedops/2.4.39 sambamba/0.8.0 subread/2.0.3 fastp/0.23.1

# fasterQ() {
#     local layout=$1
#     local sra_dir=$2
#     local sra_id=$3

#     if [ $layout = 'paired' ]; then
#         fasterq-dump $sra_dir/$sra_id.sra --split-files -O $sra_dir/

#     elif [ $layout = 'single' ]; then
#         fasterq-dump $sra_dir/$sra_id.sra -O $sra_dir/
#     fi
# }

trimming() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3

    if [ $layout = 'paired' ]; then
        fastp -w 16 \
              -i $sra_dir/raw_reads/${sra_id}_1.fastq \
              -I $sra_dir/raw_reads/${sra_id}_2.fastq \
              -o $sra_dir/trimmed_reads/${sra_id}_paired_1.fastq \
              -O $sra_dir/trimmed_reads/${sra_id}_paired_2.fastq \
              -h $sra_dir/logs/fastp.html
    elif [ $layout = 'single' ]; then
        fastp -w 16 \
              -i $sra_dir/raw_reads/${sra_id}.fastq \
              -o $sra_dir/trimmed_reads/${sra_id}_trimmed.fastq \
              -h $sra_dir/logs/fastp.html
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
             --readFilesIn $sra_dir/trimmed_reads/${sra_id}_paired_1.fastq \
                           $sra_dir/trimmed_reads/${sra_id}_paired_2.fastq \
             --outFileNamePrefix $sra_dir/logs/ \
             --outSAMtype BAM Unsorted \
             --sjdbGTFfile $host_dir/genomic.gtf \
             --outReadsUnmapped Fastx
    elif [ $layout = 'single' ]
    then
        STAR --genomeDir $genome_dir/ \
             --runThreadN 16 \
             --readFilesIn $sra_dir/trimmed_reads/${sra_id}_trimmed.fastq \
             --outFileNamePrefix $sra_dir/logs/ \
             --outSAMtype BAM Unsorted \
             --sjdbGTFfile $host_dir/genomic.gtf \
             --outReadsUnmapped Fastx
    fi
}

pathogen_alignment() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3
    local genome_dir=$4

    if [ $layout = 'paired' ]
    then
        segemehl.x -t 16 \
                   -i $genome_dir/genomic.idx \
                   -d $genome_dir/genomic.fna \
                   -q $sra_dir/trimmed_reads/${sra_id}_paired_1.fastq \
                   -p $sra_dir/trimmed_reads/${sra_id}_paired_2.fastq \
                   -o $sra_dir/sam-bam/pathogen.sam &> $sra_dir/logs/segemehl.out

    elif [ $layout = 'single' ]
    then
        segemehl.x -t 16 \
                   -i $genome_dir/genomic.idx \
                   -d $genome_dir/genomic.fna \
                   -q $sra_dir/trimmed_reads/${sra_id}_trimmed.fastq \
                   -o $sra_dir/sam-bam/pathogen.sam &> $sra_dir/logs/segemehl.out
    fi
}

infection_alignment() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3
    local genome_dir=$4

    if [ $layout = 'paired' ]
    then
        segemehl.x -t 16 \
                   -i $genome_dir/genomic.idx \
                   -d $genome_dir/genomic.fna \
                   -q $sra_dir/unmapped_reads/Unmapped.out.mate1 \
                   -p $sra_dir/unmapped_reads/Unmapped.out.mate2 \
                   -o $sra_dir/sam-bam/pathogen.sam &> $sra_dir/logs/segemehl.out

    elif [ $layout = 'single' ]
    then
        segemehl.x -t 16 \
                   -i $genome_dir/genomic.idx \
                   -d $genome_dir/genomic.fna \
                   -q $sra_dir/unmapped_reads/Unmapped.out.mate1 \
                   -o $sra_dir/sam-bam/pathogen.sam &> $sra_dir/logs/segemehl.out
    fi
}

gene_quantification() {
    local layout=$1
    local sra_dir=$2
    local sra_id=$3
    local genome_dir=$4
    local type=$5 # Indicate if it is host or pathogen

    if [ $layout = 'paired' ]
    then
        featureCounts -T 16 \
                      -p --countReadPairs \
                      -s $strandness \
                      -t gene \
                      -a $genome_dir/genomic.gtf \
                      -o $sra_dir/logs/count_matrix_${type}.txt \
                         $sra_dir/sam-bam/${type}_sorted_name.bam
    
    elif [ $layout = 'single' ]
    then
        featureCounts -T 16 \
                      -t gene \
                      -a $genome_dir/genomic.gtf \
                      -o $sra_dir/logs/count_matrix_${type}.txt \
                         $sra_dir/sam-bam/${type}_sorted_name.bam
    fi
}

## ---------------------------------- ##
## -------------  MAIN  ------------- ##
## ---------------------------------- ##
data=$1
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

for sra_dir in ~/scratch/dual_rna/experiments/$working_dir/*/*
do

    # Check if the item is a directory
    if [ -d $sra_dir ]
    then

        experiment_type="${sra_dir%/*}"
        experiment_type="${experiment_type##*/}"

        sra_id=$(basename $sra_dir)
        echo "## $sra_id"
        # Print the directory name

        # if [[ $experiment_type == 'pathogen' || $experiment_type == 'host' ]]
        # then
        #     echo "### fasterq-dump"
        #     fasterQ $layout $sra_dir $sra_id
        #     echo ""
        # fi
        
        echo "### fastp"
        trimming $layout $sra_dir $sra_id
        echo ""
        
        if [[ $experiment_type == 'host' ]]
        then
            echo "### STAR"
            host_alignment $layout $sra_dir $sra_id $host_dir
            mv $sra_dir/logs/Aligned.out.bam $sra_dir/sam-bam/host.bam
            sambamba sort -t 16 \
                            -n $sra_dir/sam-bam/host.bam \
                            -o $sra_dir/sam-bam/host_sorted_name.bam
            echo ""
        fi

        if [[ $experiment_type == 'infection' ]]
        then
            echo "### STAR"
            host_alignment $layout $sra_dir $sra_id $host_dir
            mv $sra_dir/logs/Unmapped.out.mate* $sra_dir/unmapped_reads/
            mv $sra_dir/logs/Aligned.out.bam $sra_dir/sam-bam/host.bam
            sambamba sort -t 16 \
                            -n $sra_dir/sam-bam/host.bam \
                            -o $sra_dir/sam-bam/host_sorted_name.bam
            echo ""

            echo "### segemehl"
            infection_alignment $layout $sra_dir $sra_id $pathogen_dir
            # Sort SAM file BAM file
            sambamba view -t 16 \
                            -S $sra_dir/sam-bam/pathogen.sam \
                            -f bam -o $sra_dir/sam-bam/pathogen.bam
            # Sort BAM file by name
            sambamba sort -t 16 \
                            -n $sra_dir/sam-bam/pathogen.bam \
                            -o $sra_dir/sam-bam/pathogen_sorted_name.bam
            echo ""
        fi

        if [[ $experiment_type == 'pathogen' ]]
        then
            echo "### segemehl"
            pathogen_alignment $layout $sra_dir $sra_id $pathogen_dir
            # Sort SAM file BAM file
            sambamba view -t 16 \
                            -S $sra_dir/sam-bam/pathogen.sam \
                            -f bam -o $sra_dir/sam-bam/pathogen.bam
            # Sort BAM file by name
            sambamba sort -t 16 \
                            -n $sra_dir/sam-bam/pathogen.bam \
                            -o $sra_dir/sam-bam/pathogen_sorted_name.bam
            echo ""
        fi
    fi
done
