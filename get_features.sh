#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=1GB
#SBATCH --output='../logs/MathFeature-%A.out'

# Activate the virtual environment
source /home/nima/projects/def-lpenacas/nima/newDual/env/bin/activate

python get_features.py
