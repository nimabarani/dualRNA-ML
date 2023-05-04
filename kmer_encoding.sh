#!/bin/bash
#SBATCH --time=04:30:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=50GB
#SBATCH --output='../logs/kmer-%A.out'

module load python/3.10

source ../env/bin/activate

python merege_dfs.py &> ../csvs/label.out
