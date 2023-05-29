#!/bin/bash
#SBATCH --time=15:30:00
#SBATCH --cpus-per-task=15
#SBATCH --mem-per-cpu=2GB

# Preprocessing

apptainer run -e ../tools/bioautoml.sif bash -c "python BioAutoML/BioAutoML-feature-mapping.py -fasta_train BioAutoML/data/inputs/pathogen_up.fna  BioAutoML/data/inputs/pathogen_down.fna  BioAutoML/data/inputs/pathogen_nd.fna -fasta_label_train up down nd -n_cpu 30 -output BioAutoML/data/output"
