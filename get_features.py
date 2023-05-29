import os
import subprocess
import numpy as np
import multiprocessing as mp

METHODS_DIR = '/home/nima/projects/def-lpenacas/nima/newDual/tools/MathFeature/methods'

def get_commands(label):
    input_path = f'/home/nima/projects/def-lpenacas/nima/newDual/datasets/input/pathogen_{label}_processed.fna'
    output_dir = f'/home/nima/projects/def-lpenacas/nima/newDual/datasets/output/{label}'

    commands = [
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'FourierClass_1.csv'), '-r', '1'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'FourierClass_2.csv'), '-r', '2'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'FourierClass_3.csv'), '-r', '3'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'FourierClass_4.csv'), '-r', '4'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'FourierClass_5.csv'), '-r', '5'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'FourierClass_6.csv'), '-r', '6'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'FourierClass_7.csv'), '-r', '7'],

        # Shannon - k
        ['python', os.path.join(METHODS_DIR, 'EntropyClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'EntropyClass_S5.csv'), '-k', '5', '-e', 'Shannon'],
        # Tsallis - k
        ['python', os.path.join(METHODS_DIR, 'EntropyClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'EntropyClass_T5.csv'), '-k', '5', '-e', 'Tsallis'],

        # Complex Networks (with threshold) - 12 * t
        ['python', os.path.join(METHODS_DIR, 'ComplexNetworksClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'CNC_5.csv'), '-k', '5'],
        
        # Complex Networks (without threshold - v2) - 27 * k	
        ['python', os.path.join(METHODS_DIR, 'ComplexNetworksClass-v2.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'CNCv2_5.csv'), '-k', '5'],

        
        # Other
        # # Nucleic acid composition (NAC) - 4
        # ['python', os.path.join(METHODS_DIR, 'ExtractionTechniques.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'ExT_NAC.csv'), '-t', 'NAC'],
        # # Di-nucleotide composition (DNC) - 16
        # ['python', os.path.join(METHODS_DIR, 'ExtractionTechniques.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'ExT_DNC.csv'), '-t', 'DNC'],
        # # Tri-nucleotide composition (TNC) - 64
        # ['python', os.path.join(METHODS_DIR, 'ExtractionTechniques.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'ExT_TNC.csv'), '-t', 'TNC'],
        # ORF Features or Coding Features - 10
        ['python', os.path.join(METHODS_DIR, 'CodingClass.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'CodingClass.csv')],
        # Fickett score - 2
        ['python', os.path.join(METHODS_DIR, 'FickettScore.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'FickettScore.csv'), '-seq', '1'],
        # Xmer k-Spaced Ymer Composition Frequency (kGap) - 4^X * 4^Y or 20^X * 20^Y
        ['python', os.path.join(METHODS_DIR, 'Kgap.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'Kgap_2.csv'), '-k', '2', '-bef', '1', '-aft', '2', '-seq', '1'],
        ['python', os.path.join(METHODS_DIR, 'Kgap.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'Kgap_3.csv'), '-k', '3', '-bef', '1', '-aft', '3', '-seq', '1'],
        ['python', os.path.join(METHODS_DIR, 'Kgap.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'Kgap_4.csv'), '-k', '3', '-bef', '2', '-aft', '2', '-seq', '1'],
        
        # Basic k-mer - 4^k
        ['python', os.path.join(METHODS_DIR, 'k-mers.py'), '-i', input_path, '-l', 'pathogen', '-o', os.path.join(output_dir, 'k-mers_5.csv'), '-k', '5', '-seq', '1'],
        
    ]
    return commands

def run_command(command):
    completed_process = subprocess.run(command, input=None, text=True, capture_output=True)
    return completed_process.returncode

# Get commands for the current label
labels = ['up', 'down', 'nd']
commands = []
for label in labels:
    commands.extend(get_commands(label))

# Create a process pool with 18 processes
with mp.Pool(processes=18) as pool:
    
    # Execute commands in parallel and get results
    cubes = pool.map(run_command, commands)
    
    # Print the results
    print(cubes)
