#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from mpl_toolkits import mplot3d


#%%
x=np.array([np.linspace(0,1,100),np.linspace(0,1,100)])
np.random.seed(3)
x=np.random.rand(2,100)
y=x[0,:]**2 # +x[1,:]**6    
y=np.sin(x[0,:]*2*3.14)
# np.random.shuffle(df)

ytr=y
Xtr=x
#%%
data=np.loadtxt('lorenzdata.dat')
ytr=data[0:500,3]
Xtr=data[0:500,0:3].T

ytst=data[500:1500,3]
Xtst=data[500:1500,0:3].T

#%%

nnodes=50
nvar=3

# np.random.seed(1)
def init_par():
    np.random.seed(10)
    W1=np.random.rand(nnodes,nvar)-0.5
    b1=np.random.rand(nnodes,1)-0.5
    W2=np.random.rand(1,nnodes)-0.5
    b2=np.random.rand(1,1)-0.5
    # b1=0
    # b2=0
    return W1,b1,W2,b2

def ReLU(Z):
    return np.maximum(0,Z)

def softmax(Z):
    # bb=(np.exp(Z)/np.sum(np.exp(Z)))
    # print(Z.max()) #overflows
    bb=(np.exp(Z)/sum(np.exp(Z)))
    return bb

def forward_propagation(W1,b1,W2,b2,X):
    Z1=W1.dot(X)+b1
    A1=ReLU(Z1)
    # A1=Z1
    Z2=W2.dot(A1)+b2
    # A2=softmax(Z2)
    A2=(Z2)
    return Z1,A1,Z2,A2

# def one_hot(Y):
#     ohY=np.zeros((Y.size,Y.max()+1))
#     ohY[np.arange(Y.size),Y]=1
#     ohY=ohY.T
#     return ohY

def der_ReLU(Z):
    return Z>0

def back_propagation(Z1,A1,Z2,A2,W2,X,Y):
    m=Y.size
    # print(m)
    # ohY=one_hot(Y)
    ohY=Y
    dZ2=A2-ohY
    dW2=1/m*dZ2.dot(A1.T)
    db2=1/m*np.sum(dZ2)
    dZ1=W2.T.dot(dZ2)*der_ReLU(Z1)

    dW1=1/m*dZ1.dot(X.T)
    db1=1/m*np.sum(dZ1)


    return dW1, db1, dW2, db2

def update_par(W1,b1,W2,b2,dW1,db1,dW2,db2,alpha):
    W1=W1-alpha*dW1
    b1=b1-alpha*db1
    W2=W2-alpha*dW2
    b2=b2-alpha*db2
    return W1,b1,W2,b2

def get_predictions(A2):
    # return np.argmax(A2,0)
    return A2

def get_accuracy(predictions,Y):
    # print(predictions,Y)
    # return np.sum(predictions == Y )/Y.size
    err=(predictions - Y)
    return np.sqrt(np.mean(err**2))

def gradient_descent(X,Y,iterations, alpha):
    W1,b1,W2,b2=init_par()
    resid=np.zeros((iterations,1))
    for i in range(iterations):
        Z1,A1,Z2,A2=forward_propagation(W1,b1,W2,b2,X)
        dW1, db1, dW2, db2 = back_propagation(Z1,A1,Z2,A2,W2,X,Y)
        W1,b1,W2,b2 = update_par(W1,b1,W2,b2,dW1,db1,dW2,db2,alpha)

        resid[i]=get_accuracy(get_predictions(A2),Y)

        if i%np.round(iterations/10) ==0:
            print('iteration: ', i)
            print('rms err: ',get_accuracy(get_predictions(A2),Y))
    return W1,b1,W2,b2,A2,resid

#%
W1,b1,W2,b2,A2,resid=gradient_descent(Xtr,ytr,1001,0.1)
#%%
plt.plot(Xtr[0,:],ytr.T,'x')
plt.plot(Xtr[0,:],A2.T,'x')
#%%
plt.semilogy(resid)


#%%
ax = plt.axes(projection='3d')
ax.view_init(20,150)
ax.plot3D(x[0,:].T, x[1,:].T, y.T,'.')
ax.plot3D(x[0,:].T, x[1,:].T, A2[0,:].T,'.')


# %%

from sklearn.neural_network import MLPRegressor

regr = MLPRegressor(random_state=1, max_iter=5000,
verbose=True,hidden_layer_sizes=(20,20,20),tol=1e-4,
 n_iter_no_change=100).fit(Xtr.T, ytr.T)
pr=regr.predict(Xtr.T)
print(regr.score(Xtr.T, ytr.T))
print(regr.score(Xtst.T, ytst.T))

#%%
plt.plot(ytr[0:100],'-o'); plt.plot(pr[0:100],'-x')
