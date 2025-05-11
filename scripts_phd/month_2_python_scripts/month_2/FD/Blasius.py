# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 12:25:01 2021

@author: Artur
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

n=31;
l=10;
h=l/(n+1);

J=np.zeros((2*n,2*n));
G=np.zeros(2*n);

F=np.zeros(n)+1;
H=np.zeros(n);

B=np.zeros(2*n);

for i in range(0,n):
    B[2*i]=F[i];
    B[2*i+1]=H[i];

Fa=0.0;
Ha=0.0;
Hb=1.0;

#%%
N=2*n;
for i in range(2,N-2,2):
    G[i]=2/h/h*(B[i+3]-2*B[i+1]+B[i-1])+1/2/h*(B[i+3]-B[i-1])*B[i];
    G[i+1]=1/2/h*(B[i+2]-B[i-2])-B[i+1];

i=0;
G[i]=2/h/h*(B[i+3]-2*B[i+1]+Ha)+1/2/h*(B[i+3]-Ha)*B[i];
G[i+1]=1/2/h*(B[i+2]-Fa)-B[i+1];

i=N-2;    
G[i]=2/h/h*(Hb-2*B[i+1]+B[i-1])+1/2/h*(Hb-B[i-1])*B[i];
G[i+1]=1/h*(B[i]-B[i-2])-B[i+1];

#%
for i in range(2,N-2,2):
    J[i,i]=1/2/h*(B[i+3]-B[i-1]);
    J[i,i-1]=2/h/h-1/2/h*B[i];
    J[i,i+1]=-4/h/h;
    J[i,i+3]=2/h/h+1/2/h*B[i];
    J[i+1,i+1]=-1;
    J[i+1,i+2]=1/2/h;
    J[i+1,i-2]=-1/2/h;
    

i=0;
J[i,i]=1/2/h*(B[i+3]-Ha);
#J[i,i-1]=2/h/h-1/2/h*B[i];
J[i,i+1]=-4/h/h;
J[i,i+3]=2/h/h+1/2/h*B[i];
J[i+1,i+1]=-1;
J[i+1,i+2]=1/2/h;
#J[i+1,i-2]=-1/2/h;

i=N-2;    
J[i,i]=1/2/h*(Hb-B[i-1]);
J[i,i-1]=2/h/h-1/2/h*B[i];
J[i,i+1]=-4/h/h;
#J[i,i+3]=2/h/h+1/2/h*B[i];
J[i+1,i+1]=-1;
J[i+1,i]=1/h;
J[i+1,i-2]=-1/h;    

#%%

dB=-np.linalg.inv(J)@G
B=B+dB;
print(np.linalg.norm(dB))
#%
for i in range(0,n):
    F[i]=B[2*i];
    H[i]=B[2*i+1];
    
plt.plot(np.linspace(0,l,n),H)

#%%
lam=np.linalg.eig(J);
plt.plot(np.real(lam[0])*h,np.imag(lam[0])*h,'x');
plt.xlabel('Real(Lambda)')
plt.ylabel('Imag(Lambda)')
