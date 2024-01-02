import argparse
import numpy as np
import pandas as pd
from hyperopt import Trials, STATUS_OK, tpe, fmin, STATUS_FAIL
from sklearn.model_selection import StratifiedKFold, cross_validate
import hyperparameters
from sklearn.decomposition import PCA
from mrmr_transformer import MRMRTransformer
from sklearn.preprocessing import StandardScaler

def get_data(path):
    df = pd.read_pickle(path)
    df.drop(columns="gene_id", inplace=True)
    df.columns = np.arange(df.shape[1])
    df = df.sample(frac=1, random_state=42, ignore_index=True)

    X = df.iloc[:, :-1]
    y = df.iloc[:, -1].map({"nd": 0, "up": 1, "down": 2})
    return X, y


def objective(params):

    model = params["model"](**params["params"])
    try:
        cv_results = cross_validate(model, X, y, cv=skf, scoring="f1_macro", n_jobs=-1)
        score = cv_results["test_score"].mean()
        return {"loss": -score, "status": STATUS_OK}
    except:
        return {"status": STATUS_FAIL}


def main(data_path, model, feature, nfeature):
    print(f"{model} - {feature} {nfeature}")
    global X, y, skf, best_score
    best_score = 0
    X, y = get_data(data_path)
    skf = StratifiedKFold(n_splits=10)
    trials = Trials()
    # log_file = f"{model}-{feature}.txt"
    if model == "lgbm":
        space = hyperparameters.lgbm_space
    elif model == "rf":
        space = hyperparameters.rf_space
    elif model == "xgboost":
        space = hyperparameters.xgb_space

    # Feature Selection
    if feature == "pca":
        std_scl = StandardScaler()
        X = std_scl.fit_transform(X, y)
        pca = PCA(n_components=nfeature, random_state=42)
        X = pca.fit_transform(X, y)
    elif feature == "mrmr":
        mrmr = MRMRTransformer(n_features=nfeature)
        X = mrmr.fit_transform(X, y)

    best = fmin(objective, space, algo=tpe.suggest, max_evals=1500, trials=trials)

    print(best)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Hyperparameter Tuning using HyperOpt")
    # Add arguments
    parser.add_argument("-i", "--input", required=True, help="Input Pickle")
    parser.add_argument("-m", "--model", required=True, help="Model")
    parser.add_argument("-f", "--feature", help="Feature selection method")
    parser.add_argument("-n", "--nfeature", required=True, help="Number of features")

    # Parse the arguments
    args = parser.parse_args()
    data_path = args.input
    model = args.model
    feature = args.feature
    nfeature = int(args.nfeature)
    main(data_path, model, feature, nfeature)
