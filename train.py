import argparse
import pandas as pd
import pickle
from lightgbm import LGBMClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.ensemble import (
    AdaBoostClassifier,
    GradientBoostingClassifier,
    RandomForestClassifier,
)
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.neighbors import KNeighborsClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import (
    MaxAbsScaler,
    MinMaxScaler,
    PowerTransformer,
    QuantileTransformer,
    RobustScaler,
    StandardScaler,
)
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier


def get_train_data(path):
    df: pd.DataFrame = pd.read_pickle(path)
    df = df.droplevel("nameseq").reset_index()
    df = df.sample(frac=1)

    X_train = df.iloc[:, 1:]
    X_train.columns = range(X_train.shape[1])
    y_train = df.iloc[:, 0].map({"nd": 0, "up": 1, "down": 2})
    return X_train, y_train


def get_selected_features(path, X):
    features_df = pd.read_csv(path)
    selected_features = features_df[features_df["counts"] >= 30]["feature"]
    return X[selected_features]


def train_models(X, y):
    skf = StratifiedKFold(n_splits=10)

    scorings = {
        "acc": "accuracy",
        "f1_macro": "f1_macro",
        "roc_macro": "roc_auc_ovr",
    }

    # scalers = [
    #     None,
    #     StandardScaler(),
    #     MinMaxScaler(),
    #     RobustScaler(),
    #     MaxAbsScaler(),
    #     QuantileTransformer(output_distribution='normal'),
    #     QuantileTransformer(output_distribution='uniform'),
    #     PowerTransformer(),
    # ]

    classifiers = [
        KNeighborsClassifier(),
        DecisionTreeClassifier(),
        RandomForestClassifier(),
        GradientBoostingClassifier(),
        AdaBoostClassifier(),
        LinearDiscriminantAnalysis(),
        XGBClassifier(),
        LGBMClassifier(),
    ]

    for classifier in classifiers:
        classifier_name = type(classifier).__name__
        cv_results = cross_validate(
            classifier, X, y, cv=skf, scoring=scorings, n_jobs=-1
        )
        print(f"{classifier_name}\t{cv_results['test_roc_macro'].mean():.4f}")
        print("-" * 50)


def main(data_path, features_path):
    X_train, y_train = get_train_data(data_path)
    selected_X_train = get_selected_features(features_path, X_train)
    train_models(X_train, y_train)
    print("SELECTED FEATURES")
    train_models(selected_X_train, y_train)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train different scalers and models")
    # Add arguments
    parser.add_argument("-i", "--input", required=True, help="Input Pickle")
    parser.add_argument("-f", "--features", required=True, help="Selected features")

    # Parse the arguments
    args = parser.parse_args()

    # # Access the arguments
    data_path = args.input
    features_path = args.features

    main(data_path, features_path)
