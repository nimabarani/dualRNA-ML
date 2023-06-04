import os
import numpy as np
import pandas as pd
import multiprocessing as mp

def main():
    root_directory = '/home/nima/projects/def-lpenacas/nima/newDual/datasets/output2'
    output_file = '/home/nima/projects/def-lpenacas/nima/newDual/datasets/pathogen.pkl'
    df_lists = []
    merged_df = pd.DataFrame()
    for root, dirs, files in os.walk(root_directory):
        isFirst = True
        for file in files:
            file_path = os.path.join(root, file)
            df = pd.read_csv(file_path)
            if isFirst:
                merged_df = df
                isFirst = False
                continue
            df.drop(columns=['label'], inplace=True)
            merged_df = pd.merge(df, merged_df, on='nameseq')
        if merged_df is not None:
            df_lists.append(merged_df)
    pd.concat(df_lists, ignore_index=True).to_pickle(output_file)
if __name__ == '__main__':
    main()
