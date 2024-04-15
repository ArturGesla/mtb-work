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
df=pd.read_csv("data/credit_default_sampled.csv")
print(df.head())

df.loc[df.default=='No', 'default']=0
df.loc[df.default=='Yes', 'default']=1
df.default=pd.to_numeric(df.default)
print(df.head())

model=smf.logit('default ~ income + balance',data=df)
result=model.fit()
print(result.summary())

#%%

proba=1.0/(1.0+np.exp(-result.fittedvalues))
plt.hist(proba,30)
#%%
yhat=result.predict(df[['income','balance']])
threshold=0.5
predicted_class=(yhat>threshold).astype(int)
print(predicted_class)
result.pred_table()
np.exp(result.params)
#%%
df=pd.read_csv("data/auto-mpg.csv")
print(df.head())
df['brand']=df.name.apply(lambda d:d.split(' ')[0])
print(df.head())
origin={1: 'American', 2:'European', 3:'Japanese'}
df['origin']=df.origin.apply(lambda d:origin[d])
print(df.head())
df.origin.value_counts()
pd.get_dummies(df.origin).head()
# #%%
# results=smf.ols('mpg ~ C(origin)',data=df).fit()
# results.summary()

#%
pd.get_dummies(df.origin).head()
df=df.merge(pd.get_dummies(df.origin),left_index=True,right_index=True)
print(df.head())
#%
results=smf.ols('mpg ~ Japanese + European',data=df).fit()
results.summary()
#%
df[['mpg','origin']].groupby(by='origin').mean().reset_index()

#%
df.brand.value_counts()

#%
encoder=ce.BinaryEncoder(cols=['brand'])
df=encoder.fit_transform(df)
df.head()

#%%
results=smf.ols('mpg ~ brand_0 + brand_1 + brand_2 + brand_3 + brand_4 + brand_5 ',data=df).fit()
results.summary()
#completely wrong

#%% Anova

results=smf.ols('mpg ~ C(origin)',data=df).fit()
results.summary()

aov_table=sm.stats.anova_lm(results)
print(aov_table)

from statsmodels.stats.multicomp import pairwise_tukeyhsd
from statsmodels.stats.multicomp import MultiComparison

mc=MultiComparison(df.mpg,df.origin)
results=mc.tukeyhsd()
print(results)

#%% polynomial regression
import seaborn as sns

df=pd.read_csv('data/ozone.csv')
df=df.drop(columns=['Day','Month']).rename(columns={'Solar.R':'Solar'})
# df.dropna(inplace=True)
df.head()
sns.pairplot(df)
#%%
df['Wind2']=df.Wind**2
df['Wind3']=df.Wind**3

M1=smf.ols('Ozone ~ Wind',df)
M2=smf.ols('Ozone ~ Wind + Wind2',df)
M3=smf.ols('Ozone ~ Wind+ Wind2+ Wind3',df)

M1.fit().summary()
M2.fit().summary()
M3.fit().summary()

#%% test
df=pd.read_csv("data/auto-mpg.csv")
origin={1: 'American', 2:'European', 3:'Japanese'}
df['origin']=df.origin.apply(lambda o:origin[o])
df=df.merge(pd.get_dummies(df.origin).astype(int),left_index=True,right_index= True)
# df.American=pd.to_numeric(df.American).astype(int)

print(df['American'].value_counts())

#%%
M_US0='American ~ mpg + cylinders + displacement + horsepower + weight  + acceleration'
res_US0=smf.logit(M_US0,data=df).fit()
res_US0.summary()
#%%
M_US1='American ~ mpg + cylinders + displacement  + weight  '
res_US1=smf.logit(M_US1,data=df).fit()
res_US1.summary()

yhat=res_US1.predict(df)
plt.hist(yhat)

print(res_US1.pred_table())

#%% 
df=pd.read_csv("data/titanic.csv")
df
result=smf.logit('survived ~ age + pclass + fare + sibsp + parch',data=df).fit()
result.summary()

df['fare2']=df.fare**2
df['fare3']=df.fare**3
df['fare4']=df.fare**4

M1='survived ~ fare'
M2='survived ~ fare + fare2'
M3='survived ~ fare + fare2 + fare3'
M4='survived ~ fare + fare2 + fare3 + fare4'

#%%
result=smf.logit(M3,data=df).fit()
result.summary()
print(result.pred_table())
#%%

M1='survived ~ fare + fare2 + fare3'
M2='survived ~ age + pclass + fare + sibsp + parch'
M3='survived ~ age + pclass + fare + sibsp + parch + C(sex)'

result=smf.logit(M3,data=df).fit()
result.summary()
print(result.pred_table())
tbl=result.pred_table()
acc=(tbl[0,0]+tbl[1,1])/1309
print(acc)
# %%
M1='survived ~ C(pclass)'
result=smf.logit(M1,data=df).fit()
result.summary()
