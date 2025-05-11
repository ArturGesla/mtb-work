#%% 
import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd
import scipy as scp
from scipy import stats
import statsmodels
import statsmodels.formula.api as smf
import statsmodels.api as sm


X=np.random.rand(30,3)
y=X[:,0]**2+X[:,1]**2+X[:,2]**2
plt.plot(X[:,0],y,'x')
plt.plot(X[:,1],y,'x')
plt.plot(X[:,2],y,'x')