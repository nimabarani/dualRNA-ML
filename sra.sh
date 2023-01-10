#!/bin/bash
#SBATCH --time=00:45:00
#SBATCH --output='../logs/prefetch-%A.out'

module load gcc/9.3.0 sra-toolkit/3.0.0

srid="SRR5581961"

prefetch $srid -O ~/scratch/dual_rna/raw_reads/SRR5581961
srapath $srid
