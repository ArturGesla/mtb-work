#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d


#%%
# x=np.array([np.linspace(0,1,100),np.linspace(0,1,100)]).T
np.random.seed(3); x=np.random.rand(100,2)
# y=x[:,0]**4 # +x[1,:]**6    
y=np.sin(x[:,0]*2*3.14)
X=x
# np.random.shuffle(df)

# ytr=y
# Xtr=x

#% lorenz
data=np.loadtxt('lorenzdata.dat')
np.random.shuffle(data)
y=data[0:,3]
X=data[0:,0:3]
# X=data[0:,2:3]


#%
b=y<40
y=y[b]
X=X[b,:]


ll=len(y)
ytr=y[0:np.round(0.7*ll).astype(int)]
ytst=y[np.round(0.7*ll).astype(int):]

Xtr=X[0:np.round(0.7*ll).astype(int),:]
Xtst=X[np.round(0.7*ll).astype(int):,:]





#%%
from sklearn.neural_network import MLPRegressor

regr = MLPRegressor(random_state=1, max_iter=5000,
verbose=True,hidden_layer_sizes=(30,30,30),tol=1e-6,activation='relu',
 n_iter_no_change=100).fit(Xtr, ytr)
plt.plot(regr.loss_curve_)

#%%
prtr=regr.predict(Xtr)
prtst=regr.predict(Xtst)
print("score train:",regr.score(Xtr, ytr))
print("score test :",regr.score(Xtst, ytst))

#%%
# plt.plot(ytr,'-o'); plt.plot(prtr,'-x')
# plt.plot(Xtr[:,0],ytr,'o'); plt.plot(Xtr[:,0],prtr,'x')
plt.plot(Xtr[:,2],ytr,'.'); plt.plot(Xtr[:,2],prtr,'x')
# plt.plot(Xtr,ytr,'.'); plt.plot(Xtr,prtr,'x')
# plt.plot(ytst,'-o'); plt.plot(prtst,'-x')
