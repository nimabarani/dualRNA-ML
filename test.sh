#!/bin/bash
#SBATCH --time=00:00:20

dir="/home/nima/scratch/dual_rna/raw_reads/experiment13/pathogen/SRA13234"

result="${dir%/*}"
result="${result##*/}"
# result="${result%/*/*}"


echo $result

if [[ $result == 'host' || $result == 'pathogen' ]]
then
    echo "Hi"
else
    echo "Bye"
fi