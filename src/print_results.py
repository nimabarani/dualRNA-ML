import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, make_scorer, roc_auc_score, roc_curve
from vae import VAETransformer
from sklearn.pipeline import Pipeline
from lightgbm import LGBMClassifier
from mrmr_transformer import MRMRTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import PCA
from sklearn.model_selection import StratifiedKFold, cross_validate
from sklearn.preprocessing import MinMaxScaler, StandardScaler, LabelBinarizer
from xgboost import XGBClassifier
from utils import get_data


def macro_averaged_auroc(y_true, y_pred_probs):
    """
    Compute the macro-averaged AUROC.
    
    Parameters:
    - y_true: true labels, one-hot encoded.
    - y_pred_probs: predicted probabilities for each class.
    
    Returns:
    - Macro-averaged AUROC.
    """
    
    label_binarizer = LabelBinarizer()
    label_binarizer.fit(y_true)
    y_true_bin = label_binarizer.transform(y_true)
    
    n_classes = y_true_bin.shape[1]
    auroc_per_class = []
    
    for i in range(n_classes):
        class_score = roc_auc_score(y_true_bin[:, i], y_pred_probs[:, i])
        auroc_per_class.append(class_score)
    
    return np.mean(auroc_per_class)


def macro_averaged_auprc(y_true, y_pred_probs):
    """
    Compute the macro-averaged AUPRC.

    Parameters:
    - y_true: true labels, one-hot encoded.
    - y_pred_probs: predicted probabilities for each class.

    Returns:
    - Macro-averaged AUPRC.
    """

    y_true_bin = label_binarize.transform(y_true)

    n_classes = y_true_bin.shape[1]
    auprc_per_class = []

    for i in range(n_classes):
        class_score = average_precision_score(y_true_bin[:, i], y_pred_probs[:, i])
        auprc_per_class.append(class_score)

    return np.mean(auprc_per_class)


def print_results(results):
    print(f"Accuracy\tmean: {results['test_acc'].mean():.4f}\tStd: {results['test_acc'].std():.4f}")
    print(f"Macro F1-score\tmean: {results['test_f1_macro'].mean():.4f}\tStd: {results['test_f1_macro'].std():.4f}")
    # print(f"Macro AUPRC\tmean: {results['test_macro_auprc'].mean():.4f}\tStd: {results['test_macro_auprc'].std():.4f}")
    print(f"Macro AUPRC\tmean: {results['test_macro_auroc'].mean():.4f}\tStd: {results['test_macro_auroc'].std():.4f}")
    print("-" * 30)


def main(data_path):
    X_train, y_train = get_data(data_path)
    global label_binarize
    label_binarize = LabelBinarizer().fit(y_train)

    skf = StratifiedKFold(n_splits=10)
    # macro_auprc_scorer = make_scorer(macro_averaged_auprc, needs_proba=True)
    macro_auroc_scorer = make_scorer(macro_averaged_auroc, needs_proba=True)
    
    scoring_metrics = {
        "acc": "accuracy",
        "f1_macro": "f1_macro",
    #    "macro_auprc": macro_auprc_scorer,
        "macro_auroc": macro_auroc_scorer,
    }

    xgb_clf = XGBClassifier(random_state=42, n_jobs=-1)
    rf_clf = RandomForestClassifier(random_state=42, n_jobs=-1)
    lgbm_clf = LGBMClassifier(random_state=42, n_jobs=-1)
    min_max_scl = MinMaxScaler()
    std_scl = StandardScaler()
    models = [xgb_clf, rf_clf, lgbm_clf]

    # for model in models:
    #     for n_feature in [16, 32, 64, 128, 256]:
    #         vae_tf = VAETransformer(X_train.shape[1], latent_dim=n_feature)
    #         vae_pipeline = Pipeline(
    #             [("scl", min_max_scl), ("vae", vae_tf), ("clf", model)]
    #         )

    #         vae_results = cross_validate(
    #             estimator=vae_pipeline,
    #             X=X_train,
    #             y=y_train,
    #             cv=skf,
    #             scoring=scoring_metrics,
    #             n_jobs=-1,
    #         )

    #         print(f"{type(model).__name__} - {n_feature}:")
    #         print_results(vae_results)

    # for model in models:
    #     for n_feature in [10, 20, 30, 40, 50]:
    #         pca = PCA(n_components=n_feature, random_state=42)

    #         pca_pipeline = Pipeline([("scl", std_scl), ("pca", pca), ("clf", model)])

    #         pca_results = cross_validate(
    #             estimator=pca_pipeline,
    #             X=X_train,
    #             y=y_train,
    #             cv=skf,
    #             scoring=scoring_metrics,
    #             n_jobs=-1,
    #         )

    #         print(f"{type(model).__name__} - {n_feature}:")
    #         print_results(pca_results)

    for model in models:
        for n_feature in [300, 400, 500, 600, 700]:
            mrmr = MRMRTransformer(n_features=n_feature, n_jobs=-1)
            mrmr_pipeline = Pipeline([("mrmr", mrmr), ("clf", model)])

            mrmr_results = cross_validate(
                estimator=mrmr_pipeline,
                X=X_train,
                y=y_train,
                cv=skf,
                scoring=scoring_metrics,
                n_jobs=-1,
            )

            print(f"{type(model).__name__} - {n_feature}:")
            print_results(mrmr_results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train different scalers and models")
    # Add arguments
    parser.add_argument("-i", "--input", required=True, help="Input Pickle")
    # Parse the arguments
    args = parser.parse_args()
    # Access the arguments
    data_path = args.input

    main(data_path)
