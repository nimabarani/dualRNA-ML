#!/bin/bash
#SBATCH --time=00:20:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2GB
#SBATCH --output='../logs/gff-%A.out'

module load StdEnv/2020 gffread/0.12.3

species="pathogen"
annDir="../raw_data/$species/$species"


gffread -T $annDir.gff -o $annDir.gtf
