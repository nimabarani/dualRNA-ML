import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    auc,
    average_precision_score,
    precision_recall_curve,
    roc_curve,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import label_binarize

N_CLASSES = 3
N_SPLITS = 10
RECALL_GRID = np.linspace(0, 1, 100)


def get_data(path):
    df: pd.DataFrame = pd.read_pickle(path)
    df.drop(columns="gene_id", inplace=True)
    df.columns = np.arange(df.shape[1])
    df = df.sample(frac=1, random_state=42, ignore_index=True)

    X = df.iloc[:, :-1].values
    y_true = df.iloc[:, -1].map({"nd": 0, "up": 1, "down": 2}).values
    return X, y_true


def compute_cv_prc(estimator, X, y_true):
    """Compute precision-recall curves and macro-average precision using
    cross-validation.

    Parameters
    ----------
    - classifier (estimator object): Classifier implementing 'fit' and 'predict_proba'.
    - X (array-like): Features data.
    - y_true (array-like): Target labels.

    Returns
    -------
    - all_precision (dict of lists): Precision values for each class.
    - macro_prc (list): Macro-average precision values across folds.
    """
    y_bin = label_binarize(y_true, classes=[0, 1, 2])
    cv = StratifiedKFold(n_splits=N_SPLITS)
    all_precision = {i: [] for i in range(N_CLASSES)}
    macro_prc = []

    for train, test in cv.split(X, y_true):
        estimator.fit(X[train], y_true[train])
        probas = estimator.predict_proba(X[test])

        prc_values = []

        for i in range(N_CLASSES):
            precision, recall, _ = precision_recall_curve(
                y_bin[test][:, i], probas[:, i]
            )
            all_precision[i].append(
                np.interp(RECALL_GRID, recall[::-1], precision[::-1])
            )
            prc_values.append(np.interp(RECALL_GRID, recall[::-1], precision[::-1]))

        macro_prc.append(np.mean(prc_values, axis=0))

    return all_precision, macro_prc


def plot_cv_prc(all_precision, macro_prc, name, path=None):
    """Plot the precision-recall curves for each class and
    macro-average in a cross-validation.

    Parameters
    ----------
    - all_precision (dict of lists): Precision values for each class.
    - macro_prc (list): macro-average precision values across folds.
    - name (str): Name of dataset to place in PR curve's title.
    - path (str), optional: If provided, the plot will be save in the path.

    Returns
    -------
    None
    """

    _, ax = plt.subplots(figsize=(11, 8))
    class_labels = ["ND", "UP", "DOWN"]

    mean_macro_prc = np.mean(macro_prc, axis=0)
    ax.plot(
        RECALL_GRID,
        mean_macro_prc,
        color="tab:red",
        label=f"macro-average (AUC: {auc(RECALL_GRID, mean_macro_prc):.2f})",
    )

    for i in range(N_CLASSES):
        precisions = np.array(all_precision[i])
        mean_precision = precisions.mean(axis=0)

        ax.fill_between(
            RECALL_GRID, precisions.min(axis=0), precisions.max(axis=0), alpha=0.2
        )
        ax.plot(
            RECALL_GRID,
            mean_precision,
            linestyle="--",
            linewidth=1,
            label=f"{class_labels[i]} (AUC: {auc(RECALL_GRID, mean_precision):.2f})",
        )

    ax.plot(
        [0, 1, 1],
        [1, 1, 0],
        linestyle="-.",
        color="tab:gray",
        linewidth=1,
        label=f"Perfect Classifier (AUC: {1:.2f})",
    )

    ax.plot(
        [0, 1],
        [1 / 3, 1 / 3],
        linestyle="--",
        color="k",
        linewidth=1,
        label=f"Random Classifier (AUC: {1/3:.2f})",
    )

    ax.set_title(
        f"Precision-recall curve from cross-validation of the {name} classifier."
    )
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.legend()
    plt.tight_layout()
    if path is not None:
        plt.savefig(path)
    plt.show()


def compute_pr(y_true, y_prob):
    """Compute precision-recall curves and macro-average precision.

    Parameters
    ----------
    - y_true (1-D array-like): True labels.
    - y_prob (array-like): Target probabilities/scores.

    Returns
    -------
    - precision (dict of lists): Precision values for each class.
    - recall (dict of lists): Precision values for each class.
    - average_precision (list): Macro-average precision.
    """
    y_true = label_binarize(y_true, classes=[0, 1, 2])
    class_labels = ["ND", "UP", "DOWN"]
    precision = dict()
    recall = dict()
    average_precision = dict()

    for i, class_label in enumerate(class_labels):
        precision[class_label], recall[class_label], _ = precision_recall_curve(
            y_true[:, i], y_prob[:, i]
        )
        average_precision[class_label] = average_precision_score(
            y_true[:, i], y_prob[:, i]
        )

    # Compute macro-average precision-recall curve
    precision["macro"], recall["macro"], _ = precision_recall_curve(
        y_true.ravel(), y_prob.ravel()
    )
    average_precision["macro"] = average_precision_score(
        y_true, y_prob, average="macro"
    )

    return precision, recall, average_precision


def plot_pr_curve(precision, recall, average_precision, name, path=None):
    """Plot the precision-recall curves for each class and macro-average.

    Parameters
    ----------
    - precision (dict of lists): Precision values for each class.
    - recall (dict of lists): Precision values for each class.
    - average_precision (list): Macro-average precision.
    - name (str): Name of dataset to place in PR curve's title.
    - path (str), optional: If provided, the plot will be save in the path.

    Returns
    -------
    None
    """
    class_labels = ["ND", "UP", "DOWN"]

    _, ax = plt.subplots(figsize=(10, 8))

    ax.plot(
        recall["macro"],
        precision["macro"],
        color="tab:red",
        label=f"macro-average (AUPRC = {average_precision['macro']:.2f})",
    )

    for class_label in class_labels:
        ax.plot(
            recall[class_label],
            precision[class_label],
            lw=1,
            linestyle="--",
            label=f"{class_label} (AUPRC: {average_precision[class_label]:.2f})",
        )

    ax.plot(
        [0, 1, 1],
        [1, 1, 0],
        linestyle="-.",
        color="tab:gray",
        linewidth=1,
        label=f"Perfect Classifier (AUC: {1:.2f})",
    )

    ax.plot(
        [0, 1],
        [1 / 3, 1 / 3],
        linestyle="--",
        color="k",
        linewidth=1,
        label=f"Random Classifier (AUC: {1/3:.2f})",
    )

    ax.set_title(f"Precision-Recall curve of {name} ")
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.legend()
    plt.legend(loc="upper right")
    plt.tight_layout()
    if path is not None:
        plt.savefig(path)
    plt.show()


def plot_roc_curve(y_test, y_prob, name, path):
    """
    Plot ROC curve with macro-average ROC curve for a multi-class dataset.

    Parameters
    ----------
    - y_test (array-like): Ground truth labels
    - y_prob (array-like): Predicted probabilities for each class
    - name (str): Name of dataset to place in ROC curve's title.
    - path (str), optional: If provided, the plot will be save in the path.
    """
    class_labels = ["ND", "UP", "DOWN"]
    fpr = dict()
    tpr = dict()
    roc_auc = dict()

    # Calculate ROC curve for each class
    for i in range(N_CLASSES):
        fpr[i], tpr[i], _ = roc_curve(y_test == i, y_prob[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])

    # Compute macro-average ROC curve and ROC area
    all_fpr = np.unique(np.concatenate([fpr[i] for i in range(N_CLASSES)]))

    mean_tpr = np.zeros_like(all_fpr)
    for i in range(N_CLASSES):
        mean_tpr += np.interp(all_fpr, fpr[i], tpr[i])

    mean_tpr /= N_CLASSES

    fpr["macro"] = all_fpr
    tpr["macro"] = mean_tpr
    roc_auc["macro"] = auc(fpr["macro"], tpr["macro"])

    _, ax = plt.subplots(figsize=(7, 7))
    ax.plot(
        fpr["macro"],
        tpr["macro"],
        linewidth=2,
        label=f"macro-average (AUC = {roc_auc['macro']:.2f})",
        color="tab:red",
    )

    for i in range(N_CLASSES):
        ax.plot(
            fpr[i],
            tpr[i],
            linestyle="--",
            label=f"{class_labels[i]} (AUC = {roc_auc[i]:0.2f})",
        )

    ax.plot(
        [0, 1],
        [0, 1],
        label=f"Random Classifier",
        color="black",
        linestyle="--",
        linewidth=1,
    )

    ax.plot(
        [0, 0, 1],
        [0, 1, 1],
        label=f"Perfect Classifier",
        color="tab:gray",
        scaley=False,
        linestyle="-.",
        linewidth=1,
    )

    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.set_title(f"ROC Curve for {name} Dataset")
    ax.axis("square")
    ax.legend(loc="lower right")
    plt.tight_layout()
    if path is not None:
        plt.savefig(path)
    plt.show()


def plot_cv_roc_curve(estimator, X, y_true, name, path=None):
    """Plot the precision-recall curves for each class and
    macro-average in a cross-validation.

    Parameters
    ----------
    - classifier (estimator object): Classifier implementing 'fit' and 'predict_proba'.
    - X (array-like): Features data.
    - y_true (array-like): Target labels.
    - name (str): Name of dataset to place in PR curve's title.
    - path (str), optional: If provided, the plot will be save in the path.

    Returns
    -------
    None
    """
    y_bin = label_binarize(y_true, classes=[0, 1, 2])
    cv = StratifiedKFold(n_splits=N_SPLITS)
    all_tpr = {i: [] for i in range(N_CLASSES)}
    macro_roc = []

    for train, test in cv.split(X, y_true):
        estimator.fit(X[train], y_true[train])
        probas = estimator.predict_proba(X[test])

        roc_values = []

        for i in range(N_CLASSES):
            fpr, tpr, _ = roc_curve(y_bin[test][:, i], probas[:, i])
            all_tpr[i].append(np.interp(RECALL_GRID, fpr, tpr))
            roc_values.append(np.interp(RECALL_GRID, fpr, tpr))

        macro_roc.append(np.mean(roc_values, axis=0))

    _, ax = plt.subplots(figsize=(7, 7))
    class_labels = ["ND", "UP", "DOWN"]

    FPR_GRID = np.linspace(0, 1, 100)

    mean_macro_roc = np.mean(macro_roc, axis=0)
    ax.plot(
        FPR_GRID,
        mean_macro_roc,
        color="tab:red",
        label=f"macro-average (AUC: {auc(FPR_GRID, mean_macro_roc):.2f})",
    )

    for i in range(N_CLASSES):
        tprs = np.array(all_tpr[i])
        mean_tpr = tprs.mean(axis=0)

        ax.fill_between(FPR_GRID, tprs.min(axis=0), tprs.max(axis=0), alpha=0.2)
        ax.plot(
            FPR_GRID,
            mean_tpr,
            linestyle="--",
            linewidth=1,
            label=f"{class_labels[i]} (AUC: {auc(FPR_GRID, mean_tpr):.2f})",
        )

    ax.plot(
        [0, 1],
        [0, 1],
        linestyle="--",
        color="k",
        linewidth=1,
        label="Random Classifier (AUC: 0.50)",
    )

    ax.plot(
        [0, 0, 1],
        [0, 1, 1],
        linestyle="-.",
        color="tab:gray",
        linewidth=1,
        label=f"Perfect Classifier (AUC: {1:.2f})",
    )

    ax.set_title(f"ROC curve from cross-validation of the {name} classifier.")
    ax.set_xlabel("False Positive Rate")
    ax.set_ylabel("True Positive Rate")
    ax.legend()
    plt.tight_layout()
    if path is not None:
        plt.savefig(path)
    plt.show()
