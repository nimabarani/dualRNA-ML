import argparse
from collections import Counter

import numpy as np
import pandas as pd
from reduction import svfs
from sklearn.preprocessing import StandardScaler


def get_features(train_x, train_y):
    _features = []
    fs = svfs(train_x, train_y, x_threshold=1, n_feature_threshold=100)
    reduced_data = fs.reduction()
    high_x = fs.high_rank_x()
    clean_features = high_x[reduced_data]
    dic_cls = fs.selection(high_x, reduced_data, clean_features)
    J = list(dic_cls.keys())
    selected_features = [clean_features[i] for i in J]
    _features.append(selected_features)
    return _features[0]


def main(input_path, output_path):
    NUM_ITERATION = 30
    print("START")
    df: pd.DataFrame = pd.read_pickle(input_path).reset_index(level="label")
    df = df.sample(frac=1)

    X_train = df.iloc[:, 1:]
    y_train = (
        df.iloc[:, 0].replace(["nd", "up", "down"], [0, 1, 2]).reset_index(drop=True)
    )
    print("Scaler STARTED")
    scaler = StandardScaler()
    scaled_X_train = scaler.fit_transform(X_train, y_train)
    scaled_X_train = pd.DataFrame(scaled_X_train)
    print("Scaler DONE")

    features_list = []
    results = []
    for i in range(NUM_ITERATION):
        print(f"Iteration #{i} STARTED")
        results.append(get_features(scaled_X_train, y_train))
        print(f"Iteration #{i} DONE")

    features_list = np.array(results).flatten()
    features_counts = Counter(features_list)

    features_counts_list = []
    for feature, count in features_counts.items():
        features_counts_list.append([feature, count])
        print(feature, count)

    pd.DataFrame(features_counts_list, columns=["feature", "counts"]).to_csv(
        output_path, index=False
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Feature Extraction using MathFeatures"
    )
    # Add arguments
    parser.add_argument("-i", "--input", required=True, help="Input Pickle file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")

    # Parse the arguments
    args = parser.parse_args()

    # # Access the arguments
    input_path = args.input
    output_path = args.output

    main(input_path, output_path)
