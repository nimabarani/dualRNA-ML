import argparse
import pandas as pd
from lightgbm import LGBMClassifier
from mrmr import mrmr_classif
from sklearn.model_selection import GridSearchCV, StratifiedKFold


def get_train_data(path):
    df: pd.DataFrame = pd.read_pickle(path)
    df = df.droplevel("nameseq").reset_index()
    df = df.sample(frac=1)

    X_train = df.iloc[:, 1:]
    X_train.columns = range(X_train.shape[1])
    y_train = df.iloc[:, 0].map({"nd": 0, "up": 1, "down": 2})
    return X_train, y_train


def get_selected_features(X, y):
    selected_features = mrmr_classif(X=X, y=y, K=200, n_jobs=-1)
    return X[selected_features]


def train_models(X, y):
    skf = StratifiedKFold(n_splits=10)
    parameter_grid = {
        "boosting_type": ["gbdt", "dart", "goss"],
        "num_leaves": [31, 63, 127],
        "learning_rate": [0.1, 0.01, 0.001],
        "n_estimators": [100, 200, 300],
        "min_child_samples": [20, 50, 100],
        "max_depth": [-1, 5, 10],
    }

    clf = LGBMClassifier(random_state=42)
    grid_search = GridSearchCV(
        clf, parameter_grid, cv=skf, n_jobs=-1, scoring="roc_auc_ovr"
    )
    grid_search.fit(X, y)

    print(grid_search.best_score_)
    print(grid_search.best_params_)


def main(data_path):
    X_train, y_train = get_train_data(data_path)
    selected_X_train = get_selected_features(X_train, y_train)
    train_models(selected_X_train, y_train)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train different scalers and models")
    # Add arguments
    parser.add_argument("-i", "--input", required=True, help="Input Pickle")

    # Parse the arguments
    args = parser.parse_args()
    data_path = args.input

    main(data_path)
