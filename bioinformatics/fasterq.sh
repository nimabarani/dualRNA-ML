#!/bin/bash
#SBATCH --time=01:45:00
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=12GB
#SBATCH --output='../logs/prefetch-%A.out'

module load gcc/9.3.0 sra-toolkit/3.0.0

# srid="SRR5581961"
# vdb-config --set /repository/user/main/public/root=~/scratch/
# prefetch $srid
# srapath $srid

for sra_dir in /home/nima/scratch/dual_rna/raw_reads/experiment10/infection/*
do
    sra_id="$(basename $sra_dir)"
    echo $sra_dir/$sra_id
    fasterq-dump $sra_dir/$sra_id.sra -O $sra_dir/
done
