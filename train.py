import numpy as np
import pandas as pd
import xgboost as xgb
import lightgbm as lgb
# from imblearn.combine import SMOTETomek
from sklearn.dummy import DummyClassifier
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.metrics import confusion_matrix
from sklearn.preprocessing import MinMaxScaler
from sklearn.metrics import classification_report
from sklearn.metrics import balanced_accuracy_score
from sklearn.ensemble import RandomForestClassifier
# from imblearn.under_sampling import RandomUnderSampler
from sklearn.model_selection import StratifiedKFold, cross_val_score, train_test_split
from sklearn.preprocessing import StandardScaler, MinMaxScaler, RobustScaler,\
      MaxAbsScaler, QuantileTransformer, PowerTransformer, FunctionTransformer


skf = StratifiedKFold(n_splits=10)
scalers_names = [None, 'StandardScaler', 'MinMaxScaler', 'RobustScaler', 'MaxAbsScaler', 'QuantileTransformer_normal', 'QuantileTransformer_uniform', 'PowerTransformer']
scalers = [None, StandardScaler(), MinMaxScaler(), RobustScaler(), MaxAbsScaler(), QuantileTransformer(output_distribution='normal'), QuantileTransformer(output_distribution='uniform'), PowerTransformer()]
scaler_results = [[], [], [], [], [], [], [], []]


for scaler in scalers:
    clf = make_pipeline(StandardScaler(), svm.SVC(C=1))
    cross_val_score(clf, X, y, cv=cv)
