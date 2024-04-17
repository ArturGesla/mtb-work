#%%
import numpy as np
a = np.linspace(5, 10, 11)
print(a[-3:-1])

import pandas as pd
df=pd.read_csv('class.csv')
df.head()
np.mean(df.score)
aa=df.loc[df['subject']=="Math",:]
np.mean(aa['score'])
#%%
np.mean(df.score)
#%%
aa=df.groupby('name').mean()
aa.loc[aa.score<10,:]
# %% pt 3
