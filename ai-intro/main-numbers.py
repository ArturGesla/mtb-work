#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# import sklearn
# import sklearn.datasets


# a=sklearn.datasets.load_digits()

# X=a.data
# y=a.target

# i=170
# plt.imshow(a.images[i])
# plt.title(str(y[i]))

#%%

df=pd.read_csv('mnist_train.csv')
df=np.array(df)

#%%
i=60
plt.imshow(np.reshape(df[i,1:],[28,28]))
plt.title(str(df[i,0]))
# %%
m,n=df.shape
np.random.shuffle(df)
#%%
test=df[0:1000,:].T
y=test[0,:]
X=test[1:,:]/255

train=df[1000:10000,:].T
ytr=train[0,:]
Xtr=train[1:,:]/255

#%%

def init_par():
    W1=np.random.rand(10,784)-0.5
    b1=np.random.rand(10,1)-0.5
    W2=np.random.rand(10,10)-0.5
    b2=np.random.rand(10,1)-0.5
    return W1,b1,W2,b2

def ReLU(Z):
    return np.maximum(0,Z)

def softmax(Z):
    return np.exp(Z)/np.sum(np.exp(Z))

def forward_propagation(W1,b1,W2,b2,X):
    Z1=W1.dot(X)+b1
    A1=ReLU(Z1)
    Z2=W2.dot(A1)+b2
    A2=softmax(Z2)
    return Z1,A1,Z2,A2

def one_hot(Y):
    ohY=np.zeros((Y.size,Y.max()+1))
    ohY[np.arange(Y.size),Y]=1
    ohY=ohY.T
    return ohY

def der_ReLU(Z):
    return Z>0

def back_propagation(Z1,A1,Z2,A2,W2,X,Y):
    m=Y.size
    ohY=one_hot(Y)
    dZ2=A2-ohY
    dW2=1/m*dZ2.dot(A1.T)
    db2=1/m*np.sum(dZ2)
    dZ1=W2.T.dot(dZ2)+der_ReLU(Z1)

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
    return np.argmax(A2,0)

def get_accuracy(predictions,Y):
    print(predictions,Y)
    return np.sum(predictions == Y )/Y.size

def gradient_descent(X,Y,iterations, alpha):
    W1,b1,W2,b2=init_par()
    for i in range(iterations):
        Z1,A1,Z2,A2=forward_propagation(W1,b1,W2,b2,X)
        dW1, db1, dW2, db2 = back_propagation(Z1,A1,Z2,A2,W2,X,Y)
        W1,b1,W2,b2 = update_par(W1,b1,W2,b2,dW1,db1,dW2,db2,alpha)
        if i%10 ==0:
            print('iteration: ', i)
            print('Accuracy: ',get_accuracy(get_predictions(A2),Y))
    return W1,b1,W2,b2

gradient_descent(Xtr,ytr,100,0.1)






# %%
