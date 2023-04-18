import os
import itertools
import pandas as pd
from skbio import Sequence
from skbio import DNA

# set the path to the directory containing CSV files
directory = './csvs/'

# get a list of all files in the directory
file_list = os.listdir(directory)

# filter the list to include only CSV files
csv_files = [f for f in file_list if f.endswith('.csv')]

# iterate over the list of CSV files and run your function on each one
for csv_file in csv_files:
    file_path = os.path.join(directory, csv_file)
    featureTable = pd.read_csv(file_path, sep='\t', names=['gene_id', 'sequence'])
    featureTable.set_index('gene_id', inplace=True)

    # Appending NucleotidesColumn(new features) and initialize with zeros
    kmer_products = itertools.product('ACGT', repeat=4)
    kmer_labels = [''.join(i) for i in kmer_products]
    
    sequences_df = pd.DataFrame(index=featureTable.index, columns=kmer_labels)
    
    kmers_series = featureTable['sequence'].apply(
        lambda x: Sequence(x).kmer_frequencies(k=4, overlap=True, relative=True)
    )

    sequences_df = sequences_df.apply(
        lambda row: pd.Series(kmers_series[row.name]), axis=1
    )
    sequences_df.fillna(0, inplace=True)
    
    save_path = os.path.join(directory, '4mer_'+csv_file)
    sequences_df.to_csv(save_path)