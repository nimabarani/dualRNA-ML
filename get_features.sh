#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=2GB
#SBATCH --output='../logs/MathFeature-%A.out'

# Activate the virtual environment
source /home/nima/projects/def-lpenacas/nima/newDual/env/bin/activate

python get_features.py -i $1 -o $2
