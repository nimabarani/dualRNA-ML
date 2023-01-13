#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=2GB
#SBATCH --output='../logs/fasterq-%A.out'

module load gcc/9.3.0  sra-toolkit/3.0.0

srid="SRR1714501"
fasterq-dump ../raw_reads/$srid/$srid.sra --split-files -O ../raw_reads/SRR1714501/
