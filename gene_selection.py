import argparse
import re

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# Convert FASTA to DataFrame
def fasta_to_dataframe(fasta_file):
    identifiers = []
    sequences = []
    with open(fasta_file, "r") as input_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            identifiers.append(record.id)
            sequences.append(str(record.seq))
    return pd.DataFrame({"gene_id": identifiers, "sequence": sequences})


# Convert a row to a SeqRecord
def row_to_seqrecord(row):
    sequence = Seq(row["sequence"])
    return SeqRecord(sequence, id=row["gene_id"], description="")


def main(fasta_file, gene_ids_file, output_file):
    fasta_df = fasta_to_dataframe(fasta_file)

    # Create a new column with sequence length
    fasta_df["length"] = fasta_df["sequence"].apply(len)

    # Sort by 'gene_id' and 'length' in descending order
    fasta_df.sort_values(["gene_id", "length"], ascending=[True, False], inplace=True)

    # Drop duplicates based on 'gene_id', keeping the first (longest) sequence
    fasta_df.drop_duplicates(subset="gene_id", keep="first", inplace=True)

    gene_ids_df = pd.read_csv(gene_ids_file)

    filtered_fasta_df = fasta_df[fasta_df["gene_id"].isin(gene_ids_df["gene_id"])]
    filtered_fasta_df = filtered_fasta_df[
        filtered_fasta_df["sequence"].str.contains("^[ACTGU]*$")
    ]
    
    # Apply the function to each row to get a series of SeqRecords
    seqRecord_list = filtered_fasta_df.apply(row_to_seqrecord, axis=1).tolist()

    # Write the SeqRecord objects to a FASTA file
    with open(output_file, "w") as output_handle:
        SeqIO.write(seqRecord_list, output_handle, "fasta-2line")


if __name__ == "__main__":
    # Create an argument parser
    parser = argparse.ArgumentParser(
        description="Filter genes based on the gene id file"
    )
    # Add arguments
    parser.add_argument("-i", "--input", help="Input FASTA file")
    parser.add_argument("-g", "--genes", help="Path to gene ids file")
    parser.add_argument("-o", "--output", help="Output file")

    # Parse the arguments
    args = parser.parse_args()

    # # Access the arguments
    fasta_file = args.input
    gene_ids_file = args.genes
    output_file = args.output

    main(fasta_file, gene_ids_file, output_file)
