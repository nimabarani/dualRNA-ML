import argparse
import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from mrmr_transformer import MRMRTransformer
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import MinMaxScaler
from xgboost import XGBClassifier


def get_data(path):
    df = pd.read_pickle(path)
    df.drop(columns="gene_id", inplace=True)
    df.columns = np.arange(df.shape[1])
    df = df.sample(frac=1, random_state=42, ignore_index=True)

    X = df.iloc[:, :-1]
    y = df.iloc[:, -1].map({"nd": 0, "up": 1, "down": 2})
    return X, y


def main(data_path, model, feature_selection):
    X_train, y_train = get_data(data_path)

    scl = MinMaxScaler()
    skf = StratifiedKFold(n_splits=10)

    if model == "lgbm":
        param_grid = {
            "clf__boosting_type": ["gbdt", "dart", "goss"],
            "clf__num_leaves": [31, 63, 127],
            "clf__learning_rate": [0.1, 0.01, 0.001],
            "clf__learning_rate": [0.1, 0.01, 0.001],
            "clf__n_estimators": [100, 200, 300, 400],
            "clf__min_child_samples": [20, 50, 100],
            "clf__max_depth": [-1, 5, 10],
        }
        clf = LGBMClassifier(random_state=42, n_jobs=-1)

    elif model == "rf":
        param_grid = {
            "clf__n_estimators": [100, 200, 300, 400],
            "clf__criterion": ["gini", "entropy"],
            "clf__max_depth": [None, 10, 20],
            "clf__min_samples_split": [2, 5, 10],
            "clf__min_samples_leaf": [1, 2, 4],
            "clf__max_features": ["auto", "sqrt"],
            "clf__bootstrap": [True, False],
        }
        clf = RandomForestClassifier(n_jobs=-1)

    elif model == "xgboost":
        param_grid = {
            "clf__learning_rate": [0.1, 0.01, 0.001],
            "clf__n_estimators": [100, 200, 300, 400, 500],
            "clf__max_depth": [3, 6, 9],
            "clf__min_child_weight": [1, 5, 10],
            "clf__gamma": [0, 0.1, 0.2],
            "clf__subsample": [0.8, 0.9, 1.0],
            "clf__colsample_bytree": [0.8, 0.9, 1.0],
            "clf__reg_alpha": [0, 0.1, 0.5],
            "clf__reg_lambda": [0, 0.1, 0.5],
        }
        clf = XGBClassifier(n_jobs=-1, random_state=42)

    if feature_selection == "pca":
        pca = PCA(n_components=20, random_state=42)
        pipeline = Pipeline([("scl", scl), ("fs", pca), ("clf", clf)])

    elif feature_selection == "mrmr":
        mrmr = MRMRTransformer(n_features=200)
        pipeline = Pipeline([("fs", mrmr), ("clf", clf)])

    print("Grid is Starting")
    search = GridSearchCV(pipeline, param_grid, cv=skf, scoring="f1_macro", n_jobs=-1)
    search.fit(X_train, y_train)
    print("Grid is Finished")

    print(search.best_score_)
    print(search.best_params_)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train different scalers and models")
    # Add arguments
    parser.add_argument("-i", "--input", required=True, help="Input Pickle")
    parser.add_argument("-m", "--model", required=True, help="Model")
    parser.add_argument("-f", "--feature", required=True, help="Feature Selection")

    # Parse the arguments
    args = parser.parse_args()
    data_path = args.input
    model = args.model
    feature_selection = args.feature
    main(data_path, model, feature_selection)
