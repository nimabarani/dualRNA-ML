import operator
import numpy as np
import pandas as pd
from sklearn.metrics import auc, precision_recall_curve, roc_auc_score
from sklearn.model_selection import RepeatedStratifiedKFold
import xgboost as xgb
from reduction import svfs
import multiprocessing as mp
from collections import Counter
from sklearn.dummy import DummyClassifier
# input: x_train and y_train output: selected features by svfs in a list

def get_features(train_x, train_y, th_irr=1, diff_threshold=1.7, th_red=4, k=50, alpha=50, beta=5):
    _features = []
    fs = svfs(train_x, train_y, th_irr, diff_threshold, th_red, k, alpha, beta)
    reduced_data = fs.reduction()
    high_x = fs.high_rank_x()
    clean_features = high_x[reduced_data]
    dic_cls = fs.selection(high_x, reduced_data, clean_features)
    J = list(dic_cls.keys())
    selected_features = [clean_features[i] for i in J]
    _features.append(selected_features)
    return _features[0]

def train_models(dataX, dataY, train_idx, test_dix):
    train_x, test_x = dataX.iloc[train_idx, :].copy(), dataX.iloc[test_dix, :].copy()
    train_y, test_y = dataY.iloc[train_idx].copy(), dataY.iloc[test_dix].copy()
    list_features = []
    best_roc_auc = 0
    best_idx = 0
    fs = svfs(train_x, train_y, x_threshold=1)
    reduced_data = fs.reduction()
    high_x = fs.high_rank_x()
    clean_features = high_x[reduced_data]
    dic_cls = fs.selection(high_x,reduced_data,clean_features)

    dummy_roc_list = []
    xgb_roc_list = []

    idx = 0
    for key, value in dic_cls.items():
        list_features.append(clean_features[key])
        X_train = train_x.iloc[:, list_features].copy()
        X_test = test_x.iloc[:, list_features].copy()
        
        dummy_model = DummyClassifier(strategy='stratified')
        dummy_model.fit(X_train, train_y) 
        dummy_yhat = dummy_model.predict_proba(X_test)
        dummy_roc_auc = roc_auc_score(test_y, dummy_yhat, multi_class='ovr')

        print(f'No Skill ROC AUC: {dummy_roc_auc:.4f}')

        xgb_model = xgb.XGBClassifier(objective='multi:softmax', num_class=3)
        xgb_model.fit(X_train, train_y)
        xgb_yhat = xgb_model.predict_proba(X_test)
        xgb_roc_auc = roc_auc_score(test_y, xgb_yhat, multi_class='ovr')
        print(f'XGBoost ROC AUC: {xgb_roc_auc:.4f}')

        if xgb_roc_auc > best_roc_auc:
            best_roc_auc = xgb_roc_auc
            best_idx = idx
        idx += 1
    print(f'Best ACC is: {best_roc_auc*100:.2f} for {best_idx+1} # of features')
    xgb_roc_list.append(best_roc_auc)
    return best_roc_auc

def main():
    input_file = '/home/nima/projects/def-lpenacas/nima/newDual/datasets/pathogen.csv'
    output_file = '/home/nima/projects/def-lpenacas/nima/newDual/datasets/features.out'
    df = pd.read_csv(input_file).drop(columns='nameseq')

    # Set the new column names
    df.columns = range(df.shape[1])
    acc_list = []
    features_list = []

    train_df = df.sample(frac=1, random_state=42)
    train_x = train_df.iloc[:100, :-1]
    train_y = train_df.iloc[:100, -1].replace(['up', 'down', 'nd'], [0, 1, 2])
    
    k_fold = RepeatedStratifiedKFold(n_splits=5, n_repeats=2, random_state=1)
    # print('SVFS')
    # features = get_features(train_x=train_x, train_y=train_y)
    # print('DONE')
    
    
    num_processes = 10  # Number of parallel processes to run
    num_iterations = 500  # Number of function iterations

    pool = mp.Pool(processes=num_processes)
    results = [pool.apply_async(train_models, args=(train_x, train_y, train_idx, test_dix)) for train_idx, test_dix in k_fold.split(train_x, train_y)]
    output = [p.get() for p in results]

    np.savetxt(output_file, output, fmt='%.5f')

if __name__ == "__main__":
    main()