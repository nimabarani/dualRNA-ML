import numpy as np
import pandas as pd
from reduction import svfs
# input: x_train and y_train output: selected features by svfs in a list

def get_features(train_x, train_y, th_irr=1, diff_threshold=1.7, th_red=4, k=50, alpha=50, beta=5):
    _features = []
    fs = svfs(train_x, train_y, th_irr, diff_threshold, th_red, k, alpha, beta)
    reduced_data = fs.reduction()
    high_x = fs.high_rank_x()
    clean_features = high_x[reduced_data]
    dic_cls = fs.selection(high_x, reduced_data, clean_features)
    J = list(dic_cls.keys())
    selected_features = [clean_features[i] for i in J]
    _features.append(selected_features)
    return _features[0]

def main():
    input_file = '/home/nima/projects/def-lpenacas/nima/newDual/datasets/pathogen.csv'
    output_file = '/home/nima/projects/def-lpenacas/nima/newDual/datasets/features.out'
    df = pd.read_csv(input_file).drop(columns='nameseq')

    # Set the new column names
    df.columns = range(df.shape[1])
    
    train_df = df.sample(frac=1, random_state=42)
    train_x = train_df.iloc[:, :-1]
    train_y = train_df.iloc[:, -1].replace(['up', 'down', 'nd'], [0, 1, 2])

    print('SVFS')
    features = get_features(train_x=train_x, train_y=train_y)
    print('DONE')
    
    np.savetxt(output_file, features, fmt='%d')

if __name__ == "__main__":
    main()