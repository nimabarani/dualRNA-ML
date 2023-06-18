import os
import numpy as np
import pandas as pd

def main():
    root_directory = '/home/nima/projects/def-lpenacas/nima/newDual/datasets/host/output'
    output_file = '/home/nima/projects/def-lpenacas/nima/newDual/datasets/host/host.pkl'
    df_lists = []
    merged_df = pd.DataFrame()
    for root, dirs, files in os.walk(root_directory):
        isFirst = True
        for file in files:
            file_path = os.path.join(root, file)
            df = pd.read_csv(file_path)
            df.set_index(['nameseq', 'label'], inplace=True)
            df.columns = [f'{column_name}_{file}' for column_name in df.columns]
            if isFirst:
                merged_df = df
                isFirst = False
                continue
            merged_df = pd.merge(df, merged_df, left_index=True, right_index=True)
        if not merged_df.empty:
            merged_df = merged_df.reindex(sorted(merged_df.columns), axis=1)
            df_lists.append(merged_df)
        print(merged_df.shape)
    pd.concat(df_lists).to_pickle(output_file)
if __name__ == '__main__':
    main()
