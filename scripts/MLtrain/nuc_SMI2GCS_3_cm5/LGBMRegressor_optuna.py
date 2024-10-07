### MODULES IMPORT ###
import lightgbm as lgb
import optuna

import os
import json
import numpy as np
import pandas as pd
import time
import joblib
from threadpoolctl import threadpool_limits

from sklearn import metrics
# from sklearn.model_selection import train_test_split, KFold
### END ###


def objective(trial):

    print(f'Trial {trial.number}:')

    SEARCH_PARAMS = {
        'learning_rate': trial.suggest_loguniform('learning_rate', 0.01, 0.5),
        'max_depth': trial.suggest_int('max_depth', 20, 60),
        'num_leaves': trial.suggest_int('num_leaves', 1024, 3072),
        'lambda_l1': trial.suggest_loguniform('lambda_l1', 1e-8, 10.0), #
        'lambda_l2': trial.suggest_loguniform('lambda_l2', 1e-8, 10.0), #
        'min_data_in_leaf': trial.suggest_int('min_data_in_leaf', 10, 1024),
        'feature_fraction': trial.suggest_uniform('feature_fraction', 0.4, 1.0), #
        'subsample': trial.suggest_uniform('subsample', 0.1, 1.0),
        'bagging_fraction': trial.suggest_uniform('bagging_fraction', 0.4, 1.0), #
        'bagging_freq': trial.suggest_int('bagging_freq', 1, 7), #
        'min_child_samples': trial.suggest_int('min_child_samples', 5, 100)
        }
    
    all_params = {**SEARCH_PARAMS, **STATIC_PARAMS}
    
    score = []
    mse = []
    rmse = []
    mae = []
    r2 = []
    best_score = 100000.0
    best_mse = 100000.0
    best_rmse = 100000.0
    best_mae = 100000.0
    best_r2 = 0.0
    nfolds = 5
    for i in range(1,nfolds+1):
        print(f'{i}. CV iteration', end =' ')

        train = df_train[df_train['Set'] != f'Train_fold{i}']
        train_x = np.asarray(train[atom_descriptor].tolist())
        train_y = np.asarray(train[target_type].tolist())
        
        vaild = df_train[df_train['Set'] == f'Train_fold{i}']
        valid_x = np.asarray(vaild[atom_descriptor].tolist())
        valid_y = np.asarray(vaild[target_type].tolist())

        dtrain = lgb.Dataset(train_x, label=train_y)
        dvalid = lgb.Dataset(valid_x, label=valid_y, reference=dtrain)
    
        pruning_callback = optuna.integration.LightGBMPruningCallback(trial, 'rmse', valid_name='valid')
        callbacks = [pruning_callback] if pruning_callback is not None else None

        print(all_params)
        model = lgb.train(all_params, dtrain,
                      num_boost_round=NUM_BOOST_ROUND,
                      early_stopping_rounds=EARLY_STOPPING_ROUNDS,
                      valid_sets=[dtrain, dvalid],
                      valid_names=['train', 'valid'],
                      callbacks=callbacks,
                      verbose_eval=False)
        
        preds = model.predict(valid_x, num_iteration=model.best_iteration)

        score.append(model.best_score['valid']['rmse'])
        mse.append(metrics.mean_squared_error(valid_y, preds))
        rmse.append(metrics.mean_squared_error(valid_y, preds, squared=False))
        mae.append(metrics.mean_absolute_error(valid_y, preds))
        r2.append(metrics.r2_score(valid_y, preds))
        print(f'Current Score (RMSE): {score[-1]:.4f}, Current MSE: {mse[-1]:.4f}, Current RMSE: {rmse[-1]:.4f}, Current MAE: {mae[-1]:.4f}, Current R2: {r2[-1]:.4f}')

        if score[-1] < best_score:
            best_score = score[-1]

        if mse[-1] < best_mse:
            best_mse = mse[-1]

        if rmse[-1] < best_rmse:
            best_rmse = rmse[-1]
        
        if mae[-1] < best_mae:
            best_mae = mae[-1]

        if r2[-1] > best_r2:
            best_r2 = r2[-1]
    
    print(f'\nMean Score (RMSE): {np.mean(score):.4f} +/- {np.std(score, ddof=1):.4f}  |  Best Score (RMSE): {best_score:.4f}')
    print(f'Mean MSE: {np.mean(mse):.4f} +/- {np.std(mse, ddof=1):.4f}  |  Best MSE: {best_mse:.4f}')
    print(f'Mean RMSE: {np.mean(rmse):.4f} +/- {np.std(rmse, ddof=1):.4f}  |  Best RMSE: {best_rmse:.4f}')
    print(f'Mean MAE: {np.mean(mae):.4f} +/- {np.std(mae, ddof=1):.4f}  |  Best MAE: {best_mae:.4f}')
    print(f'Mean R2: {np.mean(r2):.4f} +/- {np.std(r2, ddof=1):.4f}  |  Best R2: {best_r2:.4f}')
    print('-----------------------------------------------------------------\n')

    return np.mean(score)




if __name__ == "__main__":

    ### REGRESSOR SETUP ###
    atom_descriptor = 'nuc_SMI2GCS_3_cm5'
    target_type = 'MCA_values'
    num_cpu = 1 # <---- CHANGE NUMBER OF CPUs HERE !!!
    os.environ['MKL_NUM_THREADS'] = str(num_cpu)
    os.environ['MKL_DOMAIN_NUM_THREADS'] = 'MKL_BLAS='+str(num_cpu)
    os.environ["OMP_NUM_THREADS"] = '1'
    
    print(f'----------------REGRESSOR SETUP----------------')
    print(f'Atom Descriptor:      {atom_descriptor}')
    print(f'Target type           {target_type}')
    print(f'Number of CPUs:       {num_cpu}')
    print('\n')
    ### END ###


    ### Training/test data split ###
    # Read dataframe with precalculated GCSmmffopt descriptors and targets
    df = pd.read_csv('/groups/kemi/ree/postdoc/applications/esnuelML/data/df_nuc_x.csv.gz', index_col=0)
    df[atom_descriptor] = df.apply(lambda row: np.array(json.loads(row[atom_descriptor])), axis=1) # turn str list into real list
    df[target_type] = df.apply(lambda row: float(row[target_type]), axis=1) # turn str list into real list
    
    df_train = df[df['Set'].str.contains('Train')]
    df_test = df[df['Set'] == 'Test']

    # Clean memory
    del df

    X_test = np.asarray(df_test[atom_descriptor].tolist())
    Y_test = np.asarray(df_test[target_type].tolist())

    print('----------------DATA SET----------------')
    print(f'Training size:           {df_train.shape[0]} ({df_train.shape[0] / (df_train.shape[0] + df_test.shape[0]):.2f})')
    print(f'Held-out test size:      {df_test.shape[0]}  ({df_test.shape[0] / (df_train.shape[0] + df_test.shape[0]):.2f})')
    print('\n')
    ### END ###
    
    
    ### RUN OPTUNA FOR HYPERPARAMETER OPTIMIZATION ###
    N_TRIALS = 500
    NUM_BOOST_ROUND = 100
    EARLY_STOPPING_ROUNDS = 50

    STATIC_PARAMS = {'boosting': 'gbdt', #'gbdt', 'dart'
                    'objective': 'regression',
                    'metric': ['MAE', 'RMSE'],
                    'feature_pre_filter': False,
                    #'is_unbalance':True,
                    "verbosity": -1
                    }
    
    print('----------------OPTUNA HYPERPARAMETER OPTIMIZATION----------------')
    print('Running {} rounds of LightGBM parameter optimisation:'.format(N_TRIALS))
    
    start = time.perf_counter()
    with threadpool_limits(limits=num_cpu, user_api='openmp'):
        study = optuna.create_study(pruner=optuna.pruners.MedianPruner(n_warmup_steps=10), direction="minimize")
        study.optimize(objective, n_trials=N_TRIALS, n_jobs=num_cpu)
    finish = time.perf_counter()

    print(f'\nOptuna hyperparameter optimization finished in {round(finish-start, 2)} second(s)')
    
    joblib.dump(study, "optuna_study.pkl")
    best_params = study.best_trial.params
    all_params = {**best_params, **STATIC_PARAMS}
    print(f'Best trial parameters:\n{all_params}')
    print(f'\n')
    ### END ###


    ### TRAIN THE REGRESSION MODEL USING THE BEST OPTUNA HYPERPARAMETERS ###    
    print(f'----------------TRAINING CLASSIFIER----------------')
    score_list = []
    best_score = 100000.0
    nfolds = 5
    for i in range(1,nfolds+1):
        
        train = df_train[df_train['Set'] != f'Train_fold{i}']
        train_x = np.asarray(train[atom_descriptor].tolist())
        train_y = np.asarray(train[target_type].tolist())
        
        vaild = df_train[df_train['Set'] == f'Train_fold{i}']
        valid_x = np.asarray(vaild[atom_descriptor].tolist())
        valid_y = np.asarray(vaild[target_type].tolist())

        train_dataset = lgb.Dataset(train_x, label=train_y)
        val_dataset = lgb.Dataset(valid_x, label=valid_y, reference=train_dataset)
        
        with threadpool_limits(limits=num_cpu, user_api='openmp'):
            model = lgb.train(all_params, train_dataset, valid_sets=[train_dataset, val_dataset], verbose_eval=200,
                        num_boost_round=10000, early_stopping_rounds=250)

            score_list.append(model.best_score['valid_1']['rmse'])
            
            if score_list[-1] < best_score:
                best_score = score_list[-1]
                best_model = model

    best_model.save_model(f"final_best_model.txt", num_iteration=best_model.best_iteration)

    print(f'\n')
    print(f'----------------TRAINING RESULTS----------------')
    print(f'Mean Score (RMSE): {np.mean(score_list):.4f} +/- {np.std(score_list, ddof=1):.4f}  |  Best Score (RMSE): {best_score:.4f}')
    print(f'\n')
    ### END ###


    ### TEST THE FINAL REGRESSION MODEL ###
    final_model = lgb.Booster(model_file='final_best_model.txt')

    #-- Test data --#
    # Testing regressor - on the atomic level
    print(f'----------------Testing regressor - on the atomic level----------------')
    Y_preds = final_model.predict(X_test, num_iteration=final_model.best_iteration)
    print('Pred. MSE:', metrics.mean_squared_error(Y_test, Y_preds))
    print('Pred. RMSE:', metrics.mean_squared_error(Y_test, Y_preds, squared=False))
    print('Pred. MAE:', metrics.mean_absolute_error(Y_test, Y_preds))
    print('Pred. R2:', metrics.r2_score(Y_test, Y_preds))
    print('\n')
    ### END ###
