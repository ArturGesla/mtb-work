#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d
from sklearn.neural_network import MLPRegressor

#%%
# x=np.array([np.linspace(0,1,100),np.linspace(0,1,100)]).T
# np.random.seed(3); x=np.random.rand(100,2); y=np.sin(x[:,0]*2*3.14)
np.random.seed(3); x=np.random.rand(100,1); y=5*x; y=y-2; y=np.abs(y) #np.sin(x*2*3.14)
np.random.seed(3); x=np.random.rand(100,1); y=np.sin(2*3.14*x)+2; 
x[:,0]=np.linspace(0,1,100); y=np.sin(2*3.14*x)+2 
# y=x[:,0]**4 # +x[1,:]**6    

X=x
# np.random.shuffle(df)
plt.plot(x,y,'.')

#%% 2D
ns=10000
np.random.seed(3); X=np.random.rand(ns,2); x=X[:,0]; y=X[:,1]
z=x*(1-x)*np.cos(4*3.14*x)*np.sin(4*3.14*y**2)**2
plt.scatter(x,y,c=z)
y=z
#%

# ytr=y
# Xtr=x
#%
#%% lorenz
# data=np.loadtxt('lorenzdata.dat')
data=np.loadtxt('lorenzdata3d.dat')
# np.random.shuffle(data)
y=data[0:,3]
X=data[0:,0:3]
# X=data[0:,2:3]
# X=data[0:,2:3]
plt.plot(X,y,'.')

# # #%
b=y<40
y=y[b]
X=X[b,:]
# plt.plot(X,y,'.')





#%%

# dat=np.array([X,y])
# np.random.shuffle(dat)
# X=dat[:,0:-1]
# y=dat[:,-1]

ll=len(y)
ytr=y[0:np.round(0.7*ll).astype(int)]
ytst=y[np.round(0.7*ll).astype(int):]

Xtr=X[0:np.round(0.7*ll).astype(int),:]
Xtst=X[np.round(0.7*ll).astype(int):,:]

from sklearn.preprocessing import StandardScaler  
scaler = StandardScaler()  
# Don't cheat - fit only on training data
scaler.fit(Xtr)  
Xtr = scaler.transform(Xtr)  
# apply same transformation to test data
Xtst = scaler.transform(Xtst)  

#%%


regr = MLPRegressor(random_state=1, max_iter=50000,
verbose=True,hidden_layer_sizes=(20,20),tol=1e-4,
activation='relu',
n_iter_no_change=1000).fit(Xtr, ytr)
#%%
plt.semilogy(regr.loss_curve_)

#%%
prtr=regr.predict(Xtr)
prtst=regr.predict(Xtst)
print("score train:",regr.score(Xtr, ytr),"mean err:", np.mean(np.abs(ytr-prtr)))
print("score test :",regr.score(Xtst, ytst),"mean err:", np.mean(np.abs(ytst-prtst)))

#%%
# plt.scatter(Xtr[:,0],Xtr[:,1],c=ytr)
plt.scatter(Xtr[:,0],Xtr[:,1],c=prtr)
#%%
# plt.plot(ytr,'-o'); plt.plot(prtr,'-x')
plt.plot(Xtr[:,0],ytr,'o'); plt.plot(Xtr[:,0],prtr,'x')
# plt.plot(Xtst[:,0],ytst,'o'); plt.plot(Xtst[:,0],prtst,'x')
# plt.plot(Xtr[:,2],ytr,'.'); plt.plot(Xtr[:,2],prtr,'x')
# plt.plot(Xtr,ytr,'.'); plt.plot(Xtr,prtr,'x')
# plt.plot(ytst,'-o'); plt.plot(prtst,'-x')
#%%
plt.plot(Xtr[0:10,0],ytr[0:10],'+'); plt.plot(Xtr[0:10,0],prtr[0:10],'x')
plt.legend("train","predict")
