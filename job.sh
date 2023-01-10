#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --output='./logs/slurm-%A.out'
#SBATCH --mem-per-cpu 8G

module load  StdEnv/2020  gcc/9.3.0  sra-toolkit/2.10.8

#fasterq-dump --split-3 SRR8078974

prefetch  SRR8078974
#prefetch SRR649944  -O SRA_download/
