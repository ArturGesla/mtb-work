#%%
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
#%%

df=pd.read_csv("data/advertising.csv")
df.corr()['Sales']

model=smf.ols(formula='Sales ~ TV + Radio + Newspaper', data=df)
result=model.fit()
print(result.summary())

#%% task
df=pd.read_csv('bike+sharing+dataset/day.csv')
df=df[['season','temp','hum','windspeed','cnt']]
df.describe()

# df['cnt']=(df.cnt-np.mean(df.cnt))/np.std(df.cnt)

# df['temp']=(df.temp-np.mean(df.temp))/np.std(df.temp)
# df['hum']=(df.hum-np.mean(df.hum))/np.std(df.hum)
# df['windspeed']=(df.windspeed-np.mean(df.windspeed))/np.std(df.windspeed)
# df['season']=(df.season-np.mean(df.season))/np.std(df.season)

df.describe()

#%%
plt.clf()
plt.scatter(df['temp'],df['cnt'])
# plt.scatter(df['hum'],df['cnt'])

# plt.scatter(df['windspeed'],df['cnt'])

#%% Q1

model=smf.ols(formula='cnt ~ temp', data=df)
result=model.fit()
print(result.summary())



model=smf.ols(formula='cnt ~ hum', data=df)
result=model.fit()
print(result.summary())

model=smf.ols(formula='cnt ~ windspeed', data=df)
result=model.fit()
print(result.summary())

#%% Q2,3

seasons = {1:'spring', 2:'summer', 3:'fall', 4:'winter'}


for season in range(1,5):

    print("--"* 20)

    print("season {}".format(seasons[season]))

    for variable in ['temp', 'hum', 'windspeed']:
    # for variable in ['hum']:

        formula = "cnt ~ {}".format(variable)

        res = smf.ols(formula, data = df[df.season == season]).fit()

        # print("- R^2 for {}: {:.2f}".format(variable, res.rsquared))
        # print(res.summary())
        print("- pval")
        print(res.pvalues)


#%% Q4

model=smf.ols(formula='cnt ~ temp', data=df)
result=model.fit()
print(result.summary())

model=smf.ols(formula='cnt ~ temp + windspeed + hum', data=df)
result=model.fit()
print(result.summary())

#%% Q5

model=smf.ols(formula='cnt ~ temp + windspeed + hum', data=df)
result=model.fit()
print(result.summary())

#%% Q6

model=smf.ols(formula='cnt ~ temp + windspeed + hum + season', data=df)
result=model.fit()
print(result.summary())

#%% Q5

model=smf.ols(formula='cnt ~ temp + windspeed + hum', data=df)
result=model.fit()
print(result.summary())


rs=(result.resid-np.mean(result.resid))/np.std(result.resid)
scp.stats.kstest(rs,'norm')

#%% Q10
df['hum2']=df.hum**2
df['temp2']=df.temp**2
df['windspeed2']=df.windspeed**2

model=smf.ols(formula='cnt ~ temp + windspeed + hum', data=df)
result=model.fit()
print(result.summary())


# model=smf.ols(formula='cnt ~ temp + temp2 + windspeed + hum', data=df)
# result=model.fit()
# print(result.summary())

# model=smf.ols(formula='cnt ~ temp + windspeed + hum + hum2', data=df)
# result=model.fit()
# print(result.summary())

model=smf.ols(formula='cnt ~ temp + windspeed + windspeed2 + hum', data=df)
result=model.fit()
print(result.summary())



