#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --output='../logs/bed-%A.out'

module load bedops/2.4.39

annDir="../raw_data/host/genomic"

gtf2bed < $annDir.gtf > $annDir.bed 
