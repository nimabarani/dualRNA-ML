import pandas as pd

host_df = pd.read_csv('../csvs/train/host.csv')
host_df.set_index(keys='gene_id', inplace=True)

pathogen_df = pd.read_csv('../csvs/train/pathogen.csv')
pathogen_df.set_index(keys='gene_id', inplace=True)

train_df = pd.merge(
    host_df[host_df['label']=='up'].assign(key=1),
    pathogen_df[pathogen_df['label']=='up'].assign(key=1),
    on='key',
    how='outer',
).drop('key', axis=1)

train_df.to_csv('../csvs/train_up.csv', index=False)

train_df = pd.merge(
    host_df[host_df['label']=='down'].assign(key=1),
    pathogen_df[pathogen_df['label']=='down'].assign(key=1),
    on='key',
    how='outer',
).drop('key', axis=1)

train_df.to_csv('../csvs/train_down.csv', index=False)

train_df = pd.merge(
    host_df[host_df['label']=='up'].assign(key=1),
    pathogen_df[pathogen_df['label']=='down'].assign(key=1),
    on='key',
    how='outer',
).drop('key', axis=1)

train_df.to_csv('../csvs/train_ud.csv', index=False)

train_df = pd.merge(
    host_df[host_df['label']=='down'].assign(key=1),
    pathogen_df[pathogen_df['label']=='up'].assign(key=1),
    on='key',
    how='outer',
).drop('key', axis=1)

train_df.to_csv('../csvs/train_du.csv', index=False)

train_df = pd.merge(
    host_df[host_df['label']=='nd'].assign(key=1),
    pathogen_df[pathogen_df['label']=='nd'].assign(key=1),
    on='key',
    how='outer',
).drop('key', axis=1)

train_df.to_csv('../csvs/train_nd.csv', index=False)