#%%
from unittest import result
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
from scipy import stats
import statsmodels
#%
import statsmodels.formula.api as smf
import statsmodels.api as sm
import category_encoders as ce
#%%
# 
df=pd.read_csv("data/ozone.csv")
df=df.dropna().rename(columns = {'Solar.R':'Solar'})
df['Wind2']=np.square(df['Wind'])
df['Temp2']=np.square(df['Temp'])
print(df.head())

def RMSE(resid):
    return np.sqrt(np.square(resid).sum())/len(resid)

formulas = ['Ozone ~ Temp',
'Ozone ~ Temp + Temp2',
'Ozone ~ Wind',
'Ozone ~ Wind + Wind2',
'Ozone ~ Temp + Wind + Solar',
'Ozone ~ Temp + Temp2 + Wind + Wind2 + Solar'
]

scores = []
for formula in formulas:
    results = smf.ols(formula, df).fit()
    scores.append( { 'model': formula,
        'RMSE':RMSE(results.resid),
        'R-squared': results.rsquared} 
    )

scores = pd.DataFrame(scores)
print(scores.sort_values(by = 'RMSE').reset_index(drop = True))

#%%
# take a random sample of 70% of the dataframe
train_index = df.sample(frac = 0.7).index
# create the train and test subset
train = df.loc[df.index.isin(train_index)]
test = df.loc[~df.index.isin(train_index)]

#%
scores = []
for formula in formulas:
    results = smf.ols(formula, train).fit()
    yhat = results.predict(test)
    resid_test = yhat - test.Ozone
    scores.append( { 'model': formula,
        'RMSE_test':RMSE(resid_test),
        'RMSE_train':RMSE(results.resid),
        'R-squared': results.rsquared} 
    )

scores = pd.DataFrame(scores)
print(scores.sort_values(by = 'RMSE_test').reset_index(drop = True))
#%%
df['Wind3'] = df['Wind']**3
df['Wind4'] = df['Wind']**4
df['Temp3'] = df['Temp']**3
df['Temp4'] = df['Temp']**4
df['Solar2'] = df['Solar']**2
df['Solar3'] = df['Solar']**3
df['Solar4'] = df['Solar']**4

formulas = [
'Ozone ~ Solar',
'Ozone ~ Solar + Solar2',
'Ozone ~ Wind',
'Ozone ~ Wind + Wind2',
'Ozone ~ Temp',
'Ozone ~ Temp + Temp2',
'Ozone ~ Temp + Wind + Solar',
'Ozone ~ Temp + Temp2 + Wind + Wind2 + Solar + Solar2',
'Ozone ~ Temp + Temp2 + Temp3 + Wind + Wind2 + Wind3 + Solar + Solar2 + Solar3',
'Ozone ~ Temp + Temp2 + Temp3 + Temp4 + Wind + Wind2 + Wind3 + Wind4 + Solar + Solar2 + Solar3 + Solar4',
'Ozone ~ Temp + Temp2 + Temp3 + Temp4 + Wind + Wind2 + Wind3 + Wind4 + Solar + Solar2 + Solar3 + Solar4 + Temp * Wind + Temp * Solar + Wind * Solar',
]


scores = []
for formula in formulas:
    results = smf.ols(formula, train).fit()
    yhat = results.predict(test)
    resid_test = yhat - test.Ozone
    scores.append( { 'model': formula,
        'RMSE_test':RMSE(resid_test),
        'RMSE_train':RMSE(results.resid),
        'R-squared': results.rsquared} 
    )

scores = pd.DataFrame(scores)

plt.plot(range(0,len(scores.RMSE_test)),(scores.RMSE_test))
plt.plot(range(0,len(scores.RMSE_test)),(scores.RMSE_train))

#%%
K = 4
np.random.seed(200)
df = df.sample(frac = 1)
indexes = np.array_split(list(df.index),K)

scores = []

for formula in formulas:
    for i in range(K):
        test_index  = indexes[i]
        train_index = [idx for idx in list(df.index) if idx not in test_index]
        train = df.loc[df.index.isin(train_index)]
        test  = df.loc[~df.index.isin(train_index)]
    
        results = smf.ols(formula, train).fit()
        yhat    = results.predict(test)
        resid_test = yhat - test.Ozone
    
        scores.append( {
            'model': formula,
            'RMSE_test':RMSE(resid_test),
            'RMSE_train':RMSE(results.resid)
        })

scores = pd.DataFrame(scores)
scores = scores.groupby(by = 'model').mean().reset_index()

plt.plot(range(0,len(scores.RMSE_test)),(scores.RMSE_test))
plt.plot(range(0,len(scores.RMSE_test)),(scores.RMSE_train))

#%%
import statsmodels.formula.api as smf
import pandas as pd

df = pd.read_csv('data/credit_default_sampled.csv')
model = smf.logit('default ~ income + balance + student', data = df)
results = model.fit()
print(results.pred_table())

from sklearn.metrics import roc_curve, auc
yhat= results.predict()
false_positive_rate, true_positive_rate, thresholds = roc_curve(df['default'], yhat)

#%% test

import matplotlib.pyplot as plt

import seaborn as sns

import statsmodels.formula.api as smf

import statsmodels.api as sm

import numpy as np

import pandas as pd


from sklearn.metrics import roc_auc_score


df = pd.read_csv('data/titanic.csv')
M = 'survived ~ pclass + sex + sibsp + parch + fare + C(embarked) '

columns = ['survived', 'pclass','sex','sibsp','parch','fare','embarked']

df = df.dropna(subset = columns)

print("mean fare: {:.2f}".format(df.fare.mean()))

print(df[df.fare > 512])
#%%
seeds = [1,8,17]


for seed in seeds:

    # create the train and test subset

    np.random.seed(seed)

    train_index = df.sample(frac = 0.7).index

    train = df.loc[df.index.isin(train_index)]

    test = df.loc[~df.index.isin(train_index)]

    # train the model

    results = smf.logit(M, train).fit()

    yhat    = results.predict(test)

    # AUC score

    auc = roc_auc_score(test['survived'], yhat)
    print(results.pred_table())

    print(" seed: {} AUC: {:.2f}  ".format(seed, auc))
#%%
    df[df.age.isna()].shape
#%%

df = pd.read_csv('data/credit_default.csv')

# Convert Yes No to 1/0

df.loc[df['default'] == 'No', 'default'] = 0

df.loc[df['default'] == 'Yes', 'default'] = 1
df.default=df.default.astype(int)

# Fit a logistic regression model

results = smf.logit('default ~ student + balance + income', df).fit()
yhat = results.predict(df)

from sklearn.metrics import confusion_matrix

confusion_matrix(df['default'], yhat > 0.5)
