from hyperopt import hp
from hyperopt.pyll import scope
import numpy as np
from xgboost import XGBClassifier
from mrmr_transformer import MRMRTransformer
from sklearn.ensemble import RandomForestClassifier
from lightgbm import LGBMClassifier
from sklearn.preprocessing import MinMaxScaler, RobustScaler, StandardScaler
from sklearn.decomposition import PCA


xgb_space = hp.choice('classifier', [{
    'model': XGBClassifier,
    'params': {
        "learning_rate": hp.loguniform("learning_rate", -5, 0),
        "n_estimators": hp.choice("n_estimators", np.arange(100, 1000, 100)),
        "max_depth": hp.choice("max_depth", range(3, 15)),
        "min_child_weight": hp.quniform("min_child_weight", 1, 10, 1),
        "subsample": hp.uniform("subsample", 0.5, 1),
        "colsample_bytree": hp.uniform("colsample_bytree", 0.5, 1),
        "gamma": hp.uniform("gamma", 0, 0.5),
        "reg_alpha": hp.loguniform("reg_alpha", -5, 2),
        "reg_lambda": hp.loguniform("reg_lambda", -5, 2),
        "objective": "multi:softmax",
        "num_class": 3,
        "n_jobs": -1,
        "random_state": 42,
        "tree_method": "gpu_hist"
    }
}])



lgbm_space = hp.choice('classifier', [{
    'model': LGBMClassifier,
    'params': {
        "boosting_type": hp.choice("boosting_type", ["gbdt", "dart"]),
        "learning_rate": hp.loguniform("learning_rate", -5, 0),
        "n_estimators": hp.choice("n_estimators", np.arange(100, 1000, 100)),
        "max_depth": hp.choice("max_depth", range(-1, 15)),  # -1 means no limit
        "num_leaves": hp.choice("num_leaves", 2 ** np.arange(1, 8)),
        "min_child_samples": hp.choice("min_child_samples", range(20, 501)),
        "subsample": hp.uniform("subsample", 0.5, 1),
        "colsample_bytree": hp.uniform("colsample_bytree", 0.5, 1),
        "reg_alpha": hp.loguniform("reg_alpha", -5, 2),
        "reg_lambda": hp.loguniform("reg_lambda", -5, 2),
        "objective": "multiclass",
        "num_class": 3,
        "n_jobs": -1,
        "random_state": 42,
    }
}])


rf_space = hp.choice('classifier', [{
    'model': RandomForestClassifier,
    'params': {
        "n_estimators": hp.choice("n_estimators", np.arange(100, 1000, 100)),
        "max_depth": hp.choice("max_depth", range(3, 15)),
        "min_samples_split": hp.choice("min_samples_split", range(2, 21)),
        "min_samples_leaf": hp.choice("min_samples_leaf", range(1, 21)),
        "max_features": hp.choice("max_features", ["sqrt", "log2", None]),
        "bootstrap": hp.choice("bootstrap", [True, False]),
        "class_weight": hp.choice("class_weight", ["balanced", "balanced_subsample", None]),
        "n_jobs": -1,
        "random_state": 42,
    }
}])


