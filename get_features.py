import argparse
import os
import subprocess
import numpy as np
import multiprocessing as mp
from subprocess import Popen, PIPE, STDOUT

METHODS_DIR = '/home/nima/projects/def-lpenacas/nima/newDual/tools/MathFeature/methods'

def get_commands(input_path, output_dir, label):

    commands = [
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'FourierClass_1.csv'), '-r', '1'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'FourierClass_2.csv'), '-r', '2'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'FourierClass_3.csv'), '-r', '3'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'FourierClass_4.csv'), '-r', '4'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'FourierClass_5.csv'), '-r', '5'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'FourierClass_6.csv'), '-r', '6'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'FourierClass_7.csv'), '-r', '7'],

        # Shannon - k
        ['python', os.path.join(METHODS_DIR, 'EntropyClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'EntropyClass_S5.csv'), '-k', '4', '-e', 'Shannon'],
        # Tsallis - k
        ['python', os.path.join(METHODS_DIR, 'EntropyClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'EntropyClass_T5.csv'), '-k', '4', '-e', 'Tsallis'],

        # Complex Networks (with threshold) - 12 * t
        # ['python', os.path.join(METHODS_DIR, 'ComplexNetworksClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'CNC_5.csv'), '-k', '4'],
        
        # Complex Networks (without threshold - v2) - 27 * k	
        # ['python', os.path.join(METHODS_DIR, 'ComplexNetworksClass-v2.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'CNCv2_5.csv'), '-k', '4'],

        
        # Other
        # # Nucleic acid composition (NAC) - 4
        # ['python', os.path.join(METHODS_DIR, 'ExtractionTechniques.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'ExT_NAC.csv'), '-t', 'NAC'],
        # # Di-nucleotide composition (DNC) - 16
        # ['python', os.path.join(METHODS_DIR, 'ExtractionTechniques.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'ExT_DNC.csv'), '-t', 'DNC'],
        # # Tri-nucleotide composition (TNC) - 64
        # ['python', os.path.join(METHODS_DIR, 'ExtractionTechniques.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'ExT_TNC.csv'), '-t', 'TNC'],
        # ORF Features or Coding Features - 10
        ['python', os.path.join(METHODS_DIR, 'CodingClass.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'CodingClass.csv')],
        # Fickett score - 2
        ['python', os.path.join(METHODS_DIR, 'FickettScore.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'FickettScore.csv'), '-seq', '1'],
        # Xmer k-Spaced Ymer Composition Frequency (kGap) - 4^X * 4^Y or 20^X * 20^Y
        ['python', os.path.join(METHODS_DIR, 'Kgap.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'Kgap_2.csv'), '-k', '2', '-bef', '1', '-aft', '2', '-seq', '1'],
        ['python', os.path.join(METHODS_DIR, 'Kgap.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'Kgap_3.csv'), '-k', '3', '-bef', '1', '-aft', '3', '-seq', '1'],
        ['python', os.path.join(METHODS_DIR, 'Kgap.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'Kgap_4.csv'), '-k', '3', '-bef', '2', '-aft', '2', '-seq', '1'],
        
        # Basic k-mer - 4^k
        ['python', os.path.join(METHODS_DIR, 'k-mers.py'), '-i', input_path, '-l', label, '-o', os.path.join(output_dir, 'k-mers_5.csv'), '-k', '4', '-seq', '1'],
        
    ]
    return commands

def run_command(command):
    completed_process = subprocess.run(command, input=None, text=True, capture_output=True)
    return completed_process.returncode

def run_rckmer(script_path, input_path, output_dir, label):
    proc = subprocess.Popen(
        [
            "python",
            script_path,
            "-i",
            input_path,
            "-o",
            os.path.join(output_dir, "RC-kmer.csv"),
            "-l",
            label,
            "-t",
            "rckmer",
            "-seq",
            "1",
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout, stderr = proc.communicate("4\n".encode())

    if proc.returncode != 0:
        raise Exception(f"Script returned with error: {stderr.decode()}")

    return stdout.decode()


def main(input_path, output_dir):
    rckmer_path = os.path.join(METHODS_DIR, "ExtractionTechniques.py")
    # Get commands for the current label
    labels = ['up', 'down', 'nd']
    commands = []
    for label in labels:
        _input_path = os.path.join(input_path, f'{label}_processed.fna')
        _output_dir = os.path.join(output_dir, f'{label}')
    
        try:
            print(run_rckmer(rckmer_path, _input_path, _output_dir, label))
        except Exception as e:
            print(e)


        commands.extend(get_commands(_input_path, _output_dir, label))

    # Create a process pool with 18 processes
    with mp.Pool(processes=10) as pool:
        
        # Execute commands in parallel and get results
        return_codes = pool.map(run_command, commands)
        
        # Print the results
        print(return_codes)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Feature Extraction using MathFeatures'
    )
    # Add arguments
    parser.add_argument('-i', '--input', required=True, help='FASTA files')
    parser.add_argument('-o', '--output', required=True, help='Output directory')

    # Parse the arguments
    args = parser.parse_args()

    # # Access the arguments
    input_path = args.input
    output_dir = args.output

    main(input_path, output_dir)