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


def prepare_data(features_path, data_path):
    selected_features = pd.read_csv(features_path)['feature'].to_numpy()

    df = pd.read_pickle(data_path).reset_index(level='label')
    df = df.sample(frac=1)
    X_train = df.iloc[:, 1:]
    X_train.columns = range(X_train.shape[1])
    X_train = X_train[selected_features]
    y_train = df.iloc[:, 0]
    y_train = y_train.map({'nd': 0, 'up': 1, 'down': 2})

    return X_train, y_train


def train_model(X_train, y_train):
    skf = StratifiedKFold(
        n_splits=10,
    )

    scorings = {
        'acc': 'accuracy',
        'f1_macro': 'f1_macro',
        'roc_macro': 'roc_auc_ovr',
    }

    scalers = [
        None,
        StandardScaler(),
        MinMaxScaler(),
        RobustScaler(),
        MaxAbsScaler(),
        QuantileTransformer(output_distribution='normal'),
        QuantileTransformer(output_distribution='uniform'),
        PowerTransformer(),
    ]

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
    results = {}

    for classifier in classifiers:
        classifier_name = type(classifier).__name__
        scaler_results = {}
        for scaler in scalers:
            scaler_name = type(scaler).__name__
            model = make_pipeline(scaler, classifier)
            cv_results = cross_validate(
                model, X_train, y_train, cv=skf, scoring=scorings, n_jobs=-1
            )
            print(
                f"{classifier_name}\t{scaler_name}\t{cv_results['test_roc_macro'].mean():.4f}"
            )
            scaler_results[scaler_name] = cv_results
        results[classifier_name] = scaler_results
        print('-' * 50)
    return results


def main(data_path, features_path, output_path):
    X_train, y_train = prepare_data(
        data_path=data_path, features_path=features_path
    )
    results = train_model(X_train=X_train, y_train=y_train)
    with open(output_path, 'wb') as file:
        pickle.dump(results, file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Train different scalers and models'
    )
    # Add arguments
    parser.add_argument(
        '-i', '--input', required=True, help='Input Pickle file'
    )
    parser.add_argument(
        '-f', '--features', required=True, help='Selected features file'
    )
    parser.add_argument(
        '-o', '--output', required=True, help='Output report path'
    )

    # Parse the arguments
    args = parser.parse_args()

    # # Access the arguments
    data_path = args.input
    features_path = args.features
    output_path = args.output

    main(data_path, features_path, output_path)
