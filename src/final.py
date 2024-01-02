import pandas as pd
from joblib import dump, load
from lightgbm import LGBMClassifier
from mrmr import mrmr_classif
from sklearn.metrics import roc_auc_score


def get_train_data(path):
    df: pd.DataFrame = pd.read_pickle(path)
    df = df.droplevel("nameseq").reset_index()
    df = df.sample(frac=1)

    X_train = df.iloc[:, 1:]
    # X_train.columns = range(X_train.shape[1])
    y_train = df.iloc[:, 0].map({"nd": 0, "up": 1, "down": 2})
    return X_train, y_train


data_path = "/home/nima/projects/def-lpenacas/nima/newDual/datasets/pathogen/data.pkl"
X_train, y_train = get_train_data(data_path)
selected_features = mrmr_classif(X=X_train, y=y_train, K=200, n_jobs=-1)
selected_X_train = X_train[selected_features]
clf = LGBMClassifier(
    boosting_type="gbdt",
    learning_rate=0.01,
    max_depth=-1,
    min_child_samples=100,
    n_estimators=300,
    num_leaves=127,
    n_jobs=-1,
)
clf.fit(selected_X_train, y_train)
dump(clf, "/home/nima/projects/def-lpenacas/nima/newDual/datasets/pathogen/model.joblib") 
test_data_path = (
    "/home/nima/projects/def-lpenacas/nima/newDual/datasets/pathogen/data2.pkl"
)
X_test, y_test = get_train_data(test_data_path)
selected_X_test = X_test[selected_features]
y_hat = clf.predict_proba(selected_X_test)
result = roc_auc_score(y_test, y_hat, average='micro', multi_class="ovr")
print(result)
