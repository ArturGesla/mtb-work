#%%
rs=5
nn=30
dpth=1
from turtle import color
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d
from sklearn.neural_network import MLPRegressor

#%%
# x=np.array([np.linspace(0,1,100),np.linspace(0,1,100)]).T
# np.random.seed(3); x=np.random.rand(100,2); y=np.sin(x[:,0]*2*3.14)
# np.random.seed(3); x=np.random.rand(100,1); y=5*x; y=y-2; y=np.abs(y) #np.sin(x*2*3.14)
# np.random.seed(3); x=np.random.rand(100,1); y=np.sin(2*3.14*x)+2; 
# x[:,0]=np.linspace(0,1,100); y=np.sin(2*3.14*x)+2 
# # y=x[:,0]**4 # +x[1,:]**6    

# X=x
# # np.random.shuffle(df)
# plt.plot(x,y,'.')

# #%% 2D
# ns=10000
# np.random.seed(3); X=np.random.rand(ns,2); x=X[:,0]; y=X[:,1]
# z=x*(1-x)*np.cos(4*3.14*x)*np.sin(4*3.14*y**2)**2
# plt.scatter(x,y,c=z)
# y=z
# #%

# ytr=y
# Xtr=x
#%
#%% lorenz
# data=np.loadtxt('lorenzdata.dat')
# data=np.loadtxt('lorenzdata3d.dat')
# # np.random.shuffle(data)
# y=data[0:,3]
# X=data[0:,0:3]
# # X=data[0:,2:3]
# # X=data[0:,2:3]
# plt.plot(X,y,'.')

data=np.loadtxt('lorenzdata2d.dat')
y=data[0:,3]
X=data[0:,0:2]

y=data[0:2000,3]
X=data[0:2000,0:2]

#%%
plt.scatter(X[:,0],X[:,1],c=y)
plt.plot(7.3,7.3,'r+')
plt.text(10,7.3,"C1",color='red')
plt.colorbar()
plt.xlabel('x')
plt.ylabel('y')
plt.title("Lifetime as the function of an initial condition")
plt.savefig('lt-x0.png',dpi=300)
#%%
lts=np.sort(y)[::-1]
bb2=np.arange(1,len(lts)+1)
plt.semilogy(lts,bb2/len(bb2),'-x')
plt.ylabel('P(T>t)')
plt.xlabel('t')
plt.title("Lifetime CDF")
plt.grid()
plt.savefig('lt-cdf.png')
#%%
# # #%
# b=y<40
# y=y[b]
# X=X[b,:]
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

# rs=3; nn=30; 
depth=tuple(30 for i in range(0,dpth))

regr = MLPRegressor(random_state=rs, max_iter=5000,
verbose=True,hidden_layer_sizes=depth,tol=1e-4,
activation='relu',
n_iter_no_change=1000).fit(Xtr, ytr)
#%
plt.semilogy(regr.loss_curve_)
np.savetxt("loss-rs-"+str(rs)+"-nn-"+str(nn)+"-dpth-"+str(dpth)+".dat",regr.loss_curve_)
#%%
#%%
prtr=regr.predict(Xtr)
prtst=regr.predict(Xtst)
print("score train:",regr.score(Xtr, ytr),"mean err:", np.mean(np.abs(ytr-prtr)))
print("score test :",regr.score(Xtst, ytst),"mean err:", np.mean(np.abs(ytst-prtst)))
#%%
exit()
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
ip=0
plt.plot(Xtr[0+ip:10+ip,0],ytr[0+ip:10+ip],'+')
plt.plot(Xtr[0+ip:10+ip,0],prtr[0+ip:10+ip],'x')
plt.legend("train","predict")
#%%
plt.scatter(Xtr[:,0],Xtr[:,1],c=np.abs(ytr-prtr))
plt.colorbar()

#%%

plt.semilogy(ytr,np.abs(ytr-prtr)/ytr,'x')

#%% table


prtr=regr.predict(Xtr)
prtst=regr.predict(Xtst)
print("score train:",regr.score(Xtr, ytr),"mean err:", np.mean(np.abs(ytr-prtr)))
print("score test :",regr.score(Xtst, ytst),"mean err:", np.mean(np.abs(ytst-prtst)))
print("smaples:",len(ytst)+len(ytr))
print("mean lt:",np.mean(y))
#%%

# btr=(ytr>20)*(ytr<80)
# Xtr2=Xtr[btr,:]
# ytr2=ytr[btr]

# btst=(ytst>20)*(ytst<40)
# Xtst2=Xtst[btst,:]
# ytst2=ytst[btst]

# prtr=regr.predict(Xtr2)
# prtst=regr.predict(Xtst2)
# print("score train:",regr.score(Xtr2, ytr2),"mean err:", np.mean(np.abs(ytr2-prtr)))
# print("score test :",regr.score(Xtst2, ytst2),"mean err:", np.mean(np.abs(ytst2-prtst)))
# print("smaples:",len(ytst2)+len(ytr2))
# print("mean lt:",np.mean(ytr2))


#%%

plt.scatter(Xtr[:,0],Xtr[:,1],c=ytr)
plt.colorbar()
plt.xlabel('x_norm')
plt.ylabel('y_norm')
plt.title("Lifetime of real train data | samples:"+str(len(ytr)))
plt.savefig('lt-x0-tr.png',dpi=300)
plt.clf()

plt.scatter(Xtr[:,0],Xtr[:,1],c=prtr)
plt.colorbar()
plt.xlabel('x_norm')
plt.ylabel('y_norm')
plt.title("Lifetime of predicted train data | samples:"+str(len(ytr)))
plt.savefig('lt-x0-tr-pr.png',dpi=300)
plt.clf()

plt.scatter(Xtst[:,0],Xtst[:,1],c=ytst)
plt.colorbar()
plt.xlabel('x_norm')
plt.ylabel('y_norm')
plt.title("Lifetime of predicted test data | samples:"+str(len(ytst)))
plt.savefig('lt-x0-tst.png',dpi=300)
plt.clf()

plt.scatter(Xtst[:,0],Xtst[:,1],c=prtst)
plt.colorbar()
plt.xlabel('x_norm')
plt.ylabel('y_norm')
plt.title("Lifetime of predicted train data | samples:"+str(len(ytst)))
plt.savefig('lt-x0-tst-pr.png',dpi=300)
plt.clf()