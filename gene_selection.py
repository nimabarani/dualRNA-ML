import os
import sys
import argparse
import pandas as pd


def validate_dir(directory):
    if not (os.path.isdir(directory) | os.path.exists(directory)):
        raise argparse.ArgumentTypeError(f"Invalid directory: {directory}")
    return directory


def main(gene_csv_path, gene_list_path, output_path):
    # Read the CSV file into a pandas DataFrame
    genes_df = pd.read_csv(gene_csv_path, sep="\t", names=["gene_id", "sequence"])

    # Read the CSV file into a pandas DataFrame
    gene_ids_df = pd.read_csv(gene_list_path)

    for label in gene_ids_df["label"].unique():
        selected_genes_df = genes_df[
            genes_df["gene_id"].isin(
                gene_ids_df[gene_ids_df["label"] == label]["gene_id"]
            )
        ]

        # Sort the DataFrame by gene_id and sequence length in
        # descending order
        selected_genes_df = selected_genes_df.sort_values(
            ["gene_id", "sequence"], ascending=[True, False]
        )

        # Identify duplicated gene_ids and keep only the longest
        # sequence for each
        selected_genes_df.drop_duplicates(
            subset="gene_id", keep="first", inplace=True, ignore_index=True
        )
        print(selected_genes_df.shape)
        # Print the resulting DataFrame
        selected_genes_df.to_csv(
            f"{output_path}/genes_{label}.csv",
            index=False,
            sep="\t",
            header=None,
        )


if __name__ == "__main__":
    # Create an argument parser
    parser = argparse.ArgumentParser(
        description="Drop shortest sequences in duplicated genes."
    )
    # Add arguments
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        type=validate_dir,
        help="Genes sequences in csv format",
    )
    parser.add_argument(
        "-g", "--genes", required=True, type=validate_dir, help="Path to gene ids list"
    )
    parser.add_argument(
        "-o", "--output", required=True, default="./output.csv", help="Output directory"
    )

    # Parse the arguments
    args = parser.parse_args()

    # # Access the arguments
    gene_csv_path = args.input
    gene_list_path = args.genes
    output_path = args.output

    main(gene_csv_path, gene_list_path, output_path)
