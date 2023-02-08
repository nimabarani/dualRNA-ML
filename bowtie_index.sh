#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=16GB
#SBATCH --output='../logs/bowtie_index-%A.out'

module load bowtie2/2.4.4

genomeDir="../genomes/pathogens/GCF_000010505.1"
# mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen


genomeDir="../genomes/pathogens/GCF_017498685.1"
# mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen



genomeDir="../genomes/pathogens/GCF_000022165.1"
# mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen



genomeDir="../genomes/pathogens/GCF_000210855.2"
# mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen



genomeDir="../genomes/pathogens/GCF_000007945.1"
# mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen



genomeDir="../genomes/pathogens/GCF_000195955.2"
# mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen



genomeDir="../genomes/pathogens/GCF_002814195.1"
# mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen



genomeDir="../genomes/pathogens/GCA_002271815.1"
# mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen



genomeDir="../genomes/pathogens/GCF_000629345.1"
# mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen



genomeDir="../genomes/pathogens/GCF_004327565.1"
# mkdir $genomeDir/index/
bowtie2-build $genomeDir/genomic.fna $genomeDir/index/pathogen


