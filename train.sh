#!/bin/bash
#SBATCH --time=10:30:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=3GB
#SBATCH --output='../logs/train-%A.out'

# Activate the virtual environment
source /home/nima/projects/def-lpenacas/nima/newDual/env/bin/activate

python train.py -i $1 -f $2 -o $3