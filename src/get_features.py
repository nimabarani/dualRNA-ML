import argparse
import os
import subprocess
import numpy as np
import multiprocessing as mp
from subprocess import Popen, PIPE, STDOUT

METHODS_DIR = '/home/nima/projects/def-lpenacas/nima/newDual/tools/MathFeature/methods'

def get_commands(fasta_file, output_dir):

    commands = [
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'FourierClass_1.csv'), '-r', '1'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'FourierClass_2.csv'), '-r', '2'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'FourierClass_3.csv'), '-r', '3'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'FourierClass_4.csv'), '-r', '4'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'FourierClass_5.csv'), '-r', '5'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'FourierClass_6.csv'), '-r', '6'],
        ['python', os.path.join(METHODS_DIR, 'FourierClass.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'FourierClass_7.csv'), '-r', '7'],

        # Shannon - k
        ['python', os.path.join(METHODS_DIR, 'EntropyClass.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'EntropyClass_S5.csv'), '-k', '4', '-e', 'Shannon'],
        # Tsallis - k
        ['python', os.path.join(METHODS_DIR, 'EntropyClass.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'EntropyClass_T5.csv'), '-k', '4', '-e', 'Tsallis'],

        # Complex Networks (with threshold) - 12 * t
        # ['python', os.path.join(METHODS_DIR, 'ComplexNetworksClass.py'), '-i', fasta_file, '-l', label, '-o', os.path.join(output_dir, 'CNC_5.csv'), '-k', '4'],
        
        # Complex Networks (without threshold - v2) - 27 * k	
        ['python', os.path.join(METHODS_DIR, 'ComplexNetworksClass-v2.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'CNCv2.csv'), '-k', '4'],

        
        # Other
        # ORF Features or Coding Features - 10
        ['python', os.path.join(METHODS_DIR, 'CodingClass.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'CodingClass.csv')],
        # Fickett score - 2
        ['python', os.path.join(METHODS_DIR, 'FickettScore.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'FickettScore.csv'), '-seq', '1'],
        # Xmer k-Spaced Ymer Composition Frequency (kGap) - 4^X * 4^Y or 20^X * 20^Y
        ['python', os.path.join(METHODS_DIR, 'Kgap.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'k-gap_1.csv'), '-k', '2', '-bef', '1', '-aft', '2', '-seq', '1'],
        ['python', os.path.join(METHODS_DIR, 'Kgap.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'k-gap_2.csv'), '-k', '3', '-bef', '1', '-aft', '3', '-seq', '1'],
        ['python', os.path.join(METHODS_DIR, 'Kgap.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'k-gap_3.csv'), '-k', '3', '-bef', '2', '-aft', '2', '-seq', '1'],
        
        # Basic k-mer - 4^k
        ['python', os.path.join(METHODS_DIR, 'k-mers.py'), '-i', fasta_file, '-o', os.path.join(output_dir, 'k-mers.csv'), '-k', '4', '-seq', '1'],
        
    ]
    return commands

def run_command(command):
    completed_process = subprocess.run(command, input=None, text=True, capture_output=True)
    return completed_process.returncode

def run_rckmer(script_path, fasta_file, output_dir):
    proc = subprocess.Popen(
        [
            "python",
            script_path,
            "-i",
            fasta_file,
            "-o",
            os.path.join(output_dir, "RC-kmer.csv"),
            "-t",
            "rckmer",
            "-seq",
            "1",
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    stdout, stderr = proc.communicate("4\n".encode())

    if proc.returncode != 0:
        raise Exception(f"Script returned with error: {stderr.decode()}")

    return stdout.decode()


def main(fasta_file, output_dir):
    rckmer_path = os.path.join(METHODS_DIR, "ExtractionTechniques.py")
    
    commands = get_commands(fasta_file, output_dir)
    
    try:
        print(run_rckmer(rckmer_path, fasta_file, output_dir))
    except Exception as e:
        print(e)


    with mp.Pool(processes=14) as pool:
        return_codes = pool.map(run_command, commands)
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
    fasta_file = args.input
    output_dir = args.output

    main(fasta_file, output_dir)