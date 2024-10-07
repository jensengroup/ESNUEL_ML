### MODULES IMPORT ###
import os
import json
from joblib import dump, load
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from sklearn import metrics
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
### END ###


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


    ### TRAIN THE REGRESSION MODEL USING THE BEST OPTUNA HYPERPARAMETERS ###    
    print(f'----------------TRAINING CLASSIFIER----------------')

    model = RandomForestRegressor(n_estimators=200)

    rmse = []
    mse = []
    mae = []
    r2 = []
    best_rmse = 100000.0
    best_mse = 100000.0
    best_mae = 100000.0
    best_r2 = 0.0
    nfolds = 5
    for i in range(1,nfolds+1):
        
        train = df_train[df_train['Set'] != f'Train_fold{i}']
        train_x = np.asarray(train[atom_descriptor].tolist())
        train_y = np.asarray(train[target_type].tolist())
        
        vaild = df_train[df_train['Set'] == f'Train_fold{i}']
        valid_x = np.asarray(vaild[atom_descriptor].tolist())
        valid_y = np.asarray(vaild[target_type].tolist())
        
        with threadpool_limits(limits=num_cpu, user_api='openmp'):
            model.fit(train_x, train_y)

            preds = model.predict(valid_x)

            mse.append(metrics.mean_squared_error(valid_y, preds))
            rmse.append(metrics.mean_squared_error(valid_y, preds, squared=False))
            mae.append(metrics.mean_absolute_error(valid_y, preds))
            r2.append(metrics.r2_score(valid_y, preds))
            print(f'Current RMSE: {rmse[-1]:.4f}, Current MSE: {mse[-1]:.4f}, Current MAE: {mae[-1]:.4f}, Current R2: {r2[-1]:.4f}')

            if rmse[-1] < best_rmse:
                best_rmse = rmse[-1]
                best_model = model

            if mse[-1] < best_mse:
                best_mse = mse[-1]

            if mae[-1] < best_mae:
                best_mae = mae[-1]

            if r2[-1] > best_r2:
                best_r2 = r2[-1]

    dump(best_model, 'final_best_model.joblib')

    print(f'\n')
    print(f'\nMean RMSE: {np.mean(rmse):.4f} +/- {np.std(rmse, ddof=1):.4f}  |  Best RMSE: {best_rmse:.4f}')
    print(f'Mean MSE: {np.mean(mse):.4f} +/- {np.std(mse, ddof=1):.4f}  |  Best MSE: {best_mse:.4f}')
    print(f'Mean MAE: {np.mean(mae):.4f} +/- {np.std(mae, ddof=1):.4f}  |  Best MAE: {best_mae:.4f}')
    print(f'Mean R2: {np.mean(r2):.4f} +/- {np.std(r2, ddof=1):.4f}  |  Best R2: {best_r2:.4f}')
    print('-----------------------------------------------------------------\n')
    ### END ###


    ### TEST THE FINAL REGRESSION MODEL ###
    final_model = load('final_best_model.joblib')

    # Testing regressor - on the atomic level
    print(f'----------------Testing regressor - on the atomic level----------------')
    Y_preds = final_model.predict(X_test)
    print('Pred. MSE:', metrics.mean_squared_error(Y_test, Y_preds))
    print('Pred. RMSE:', metrics.mean_squared_error(Y_test, Y_preds, squared=False))
    print('Pred. MAE:', metrics.mean_absolute_error(Y_test, Y_preds))
    print('Pred. R2:', metrics.r2_score(Y_test, Y_preds))
    print('\n')
    ### END ###
