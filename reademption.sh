#!/bin/bash
#SBATCH --time=8:30:00
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=8GB
#SBATCH --output='../logs/reademption-%A.out'


source ../env/bin/activate

reademption align -f /home/nima/projects/def-lpenacas/nima/newDual/scratch/READemption_analysis \
                  -p 10 \
                  -P \
                  -q \
                  -g