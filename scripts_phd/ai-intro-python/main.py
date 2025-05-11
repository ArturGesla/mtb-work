# 
import pandas as pd
#import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp
#%%
df=pd.read_csv('bike+sharing+dataset/day.csv')
df = df[['season', 'temp','hum','windspeed','cnt']]

#sns.pairplot(df)

#plt.plot(df['season'])
plt.clf();
#plt.scatter(df['hum'],df['temp'])
plt.scatter(df['cnt'],df['temp'])
#plt.scatter(df['windspeed'],df['hum'])

df.corr()

#%%

b=(df['season']==4); np.mean(df['cnt'][b])

#%%

b=(df['season']==3); us3=(df['cnt'][b])
b=(df['season']==2); us2=(df['cnt'][b])

scp.stats.ttest_ind(us3,us2)

#%%

df['cnt']=(df.cnt-np.mean(df.cnt))/np.std(df.cnt)
df['cnt'].describe()
scp.stats.kstest(df.cnt,'norm')

#%%
df['hum']=(df.hum-np.mean(df.hum))/np.std(df.hum)
df['hum'].describe()
scp.stats.kstest(df.hum,'norm')

#%%
df['hum']=(df.hum-np.mean(df.hum))/np.std(df.hum)
df['hum'].describe()
scp.stats.kstest(df.hum,'norm')




