import sys
from debiasm import DebiasMClassifier
import debiasm
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

np.random.seed(123)
log_DEBIAS = pd.read_csv("./csv_files/autism_log_DEBIAS-M.csv") # read in data
data = log_DEBIAS.iloc[:,3:len(log_DEBIAS)] # select the case, Study ID (0-14), and data columns

## the read count matrix and "batches" a.k.a sample_IDs
X = data.iloc[:,2:len(data)]
X = X.to_numpy()

## select the labels
y =  data["case"].to_numpy() # the 0s and 1s... make sure true (case/1) false(control/0)

# batches
# batches = data.iloc[:,1:2]
# batches = batches.to_numpy()

# batches
batches = data.iloc[:, 1]
batches  = batches.to_numpy().astype(int)


## we assume the batches are numbered ints starting at '0',
## and they are in the first column of the input X matrices

# run this one without editing
# X_with_batch = np.hstack((batches[:, np.newaxis], X))
X_with_batch = np.hstack((batches[:, np.newaxis], X))
X_with_batch[:10, :10]


## there will be 15 different validations with 15 testing data resulting in 15 different spreadsheets

# this means the training will be done without the fourth batch
i = int(sys.argv[1])

val_inds = batches==i
X_train, X_val = X_with_batch[~val_inds], X_with_batch[val_inds]
y_train, y_val = y[~val_inds], y[val_inds]


dmc = DebiasMClassifier(x_val=X_val) ## give it the held-out inputs to account for
                                    ## those domains shifts while training

dmc.fit(X_train, y_train)
print('finished training!')

print(roc_auc_score(y_val, dmc.predict_proba(X_val)[:, 1]))

X_debiassed = dmc.transform(X_with_batch)


#### write to CSV ####
id_dict = {}
Study_ID = log_DEBIAS["Study_ID"].unique()
id_array = log_DEBIAS["ID"].unique()
for id in range(0,len(Study_ID)):
    id_dict[id_array[id]] = Study_ID[id]

cols = list(data.columns[2:len(data)])

debiased_df = pd.DataFrame(data = X_debiassed,
                           columns = cols)

first_cols = log_DEBIAS.iloc[:, 0:3]

frames = [first_cols, debiased_df]

merged = pd.concat(frames, axis=1)

merged.to_csv(f"./csv_files/DEBIAS-M_runs/semi_leaky_autism/autism_debiased_lognorm_{id_dict[i]}.csv", index = False)


weights = dmc.model.batch_weights
np_weights = weights.detach().numpy()
cols = list(data.columns[2:len(data)])


# Study_ID = log_DEBIAS["Study_ID"].unique()

weights_df = pd.DataFrame(data = np_weights,
                          columns = cols,
                          index = Study_ID)

weights_df = weights_df.transpose()



weights_df.to_csv(f"./csv_files/DEBIAS-M_runs/semi_leaky_autism/autism_debias_weights_{id_dict[i]}.csv")
