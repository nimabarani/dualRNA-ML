#!/bin/bash
#SBATCH --time=05:30:00
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2GB
#SBATCH --output='../logs/SVFS-%A.out'

# Activate the virtual environment
source /home/nima/projects/def-lpenacas/nima/newDual/env/bin/activate

python ./svfs/main.py -i $1 -o $2