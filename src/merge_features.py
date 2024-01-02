import argparse
from pathlib import Path

import pandas as pd


def main(input_dir, gene_ids_file, output_path):
    csv_files = Path(input_dir)

    gene_ids_df = pd.read_csv(gene_ids_file)
    merged_df = pd.DataFrame()

    isFirst = True
    for csv_file in csv_files.iterdir():
        df = pd.read_csv(csv_file).drop(columns="label")
        df.set_index("nameseq", inplace=True)
        df.columns = [f"{column_name}_{csv_file.stem}" for column_name in df.columns]

        if isFirst:
            merged_df = df
            isFirst = False
            continue

        merged_df = pd.merge(df, merged_df, left_index=True, right_index=True)

    if not merged_df.empty:
        merged_df = merged_df.reindex(sorted(merged_df.columns), axis=1)

    merged_df = merged_df.reset_index(names="gene_id").merge(gene_ids_df, on="gene_id")
    merged_df.to_pickle(output_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge features")
    # Add arguments
    parser.add_argument("-i", "--input", required=True, help="CSV files directory")
    parser.add_argument("-g", "--genes", help="Path to gene ids file")
    parser.add_argument("-o", "--output", required=True, help="Final pickle path")

    # Parse the arguments
    args = parser.parse_args()

    # # Access the arguments
    input_dir = args.input
    gene_ids_file = args.genes
    output_path = args.output

    main(input_dir, gene_ids_file, output_path)
