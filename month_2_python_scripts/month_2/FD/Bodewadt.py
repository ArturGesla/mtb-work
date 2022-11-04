# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 12:25:01 2021

@author: Artur
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

Re=1;

n=99;
l=14;
h=l/(n+1);
fields=3;

J=np.zeros((fields*n,fields*n));
g=np.zeros(fields*n);

F=np.zeros(n);
G=np.zeros(n)+1;
H=np.zeros(n)+1;

B=np.zeros(fields*n);

for i in range(0,n):
    B[fields*i]=F[i];
    B[fields*i+1]=G[i];
    B[fields*i+2]=H[i];

Fa=0.0;
Ga=0.0;
Ha=0.0;
Fb=0.0;
Gb=1.0;
#Hb=1.0;

#%%
N=fields*n;
for i in range(3,N-3,fields):
    g[i]=1/2/h*(B[i+5]-B[i-1])+2*B[i];
    g[i+1]=1/h/h*(B[i+3]-2*B[i]+B[i-3])+B[i+1]*B[i+1]-B[i]*B[i]-1/2/h*(B[i+3]-B[i-3])*B[i+2]-1;
    g[i+2]=1/h/h*(B[i+4]-2*B[i+1]+B[i-2])-1/2/h*(B[i+4]-B[i-2])*B[i+2]-2*B[i]*B[i+1];

i=0;
g[i]=1/2/h*(B[i+5]-Ha)+2*B[i];
g[i+1]=1/h/h*(B[i+3]-2*B[i]+Fa)+B[i+1]*B[i+1]-B[i]*B[i]-1/2/h*(B[i+3]-Fa)*B[i+2]-1;
g[i+2]=1/h/h*(B[i+4]-2*B[i+1]+Ga)-1/2/h*(B[i+4]-Ga)*B[i+2]-2*B[i]*B[i+1];


i=N-3;    
g[i]=1/h*(B[i+2]-B[i-1])+2*B[i]; # Diferent discretisation!!!
g[i+1]=1/h/h*(Fb-2*B[i]+B[i-3])+B[i+1]*B[i+1]-B[i]*B[i]-1/2/h*(Fb-B[i-3])*B[i+2]-1;
g[i+2]=1/h/h*(Gb-2*B[i+1]+B[i-2])-1/2/h*(Gb-B[i-2])*B[i+2]-2*B[i]*B[i+1];


#%
for i in range(3,N-3,3):
    J[i,i]=2;
    J[i,i-1]=-1/2/h;
    J[i,i+5]=1/2/h;
    
    J[i+1,i]=-2/h/h-2*B[i];
    J[i+1,i+1]=2*B[i+1];
    J[i+1,i+2]=-1/2/h*(B[i+3]-B[i-3]);
    J[i+1,i+3]=1/h/h-1/2/h*B[i+2];
    J[i+1,i-3]=1/h/h+1/2/h*B[i+2];

    
    J[i+2,i]=-2*B[i+1];
    J[i+2,i+1]=-2*B[i]-2/h/h;
    J[i+2,i+2]=-1/2/h*(B[i+4]-B[i-2]);
    J[i+2,i+4]=1/h/h-1/2/h*B[i+2];
    J[i+2,i-2]=1/h/h+1/2/h*B[i+2];
    
    
i=0;
J[i,i]=2;
#J[i,i-1]=-1/2/h;
J[i,i+5]=1/2/h;

J[i+1,i]=-2/h/h-2*B[i];
J[i+1,i+1]=2*B[i+1];
J[i+1,i+2]=-1/2/h*(B[i+3]-Fa);
J[i+1,i+3]=1/h/h-1/2/h*B[i+2];
#J[i+1,i-3]=1/h/h+1/2/h*B[i+2];


J[i+2,i]=-2*B[i+1];
J[i+2,i+1]=-2*B[i]-2/h/h;
J[i+2,i+2]=-1/2/h*(B[i+4]-Ga);
J[i+2,i+4]=1/h/h-1/2/h*B[i+2];
#J[i+2,i-2]=1/h/h+1/2/h*B[i+2];

i=N-3;    
J[i,i]=2;
J[i,i-1]=-1/h;
J[i,i+2]=1/h;

J[i+1,i]=-2/h/h-2*B[i];
J[i+1,i+1]=2*B[i+1];
J[i+1,i+2]=-1/2/h*(Fb-B[i-3]);
#J[i+1,i+3]=1/h/h-1/2/h*B[i+2];
J[i+1,i-3]=1/h/h+1/2/h*B[i+2];


J[i+2,i]=-2*B[i+1];
J[i+2,i+1]=-2*B[i]-2/h/h;
J[i+2,i+2]=-1/2/h*(Gb-B[i-2]);
#J[i+2,i+4]=1/h/h-1/2/h*B[i+2];
J[i+2,i-2]=1/h/h+1/2/h*B[i+2];

#%

dB=-np.linalg.inv(J)@g
B=B+dB;
print(np.linalg.norm(dB))
#%%
for i in range(0,n):
    F[i]=B[fields*i];
    G[i]=B[fields*i+1];
    H[i]=B[fields*i+2];
    
plt.plot(np.linspace(0,l,n+2),np.append(np.append(Fa,F),Fb))
plt.plot(np.linspace(0,l,n+2),np.append(np.append(Ga,G),Gb))
plt.plot(np.linspace(0,l,n+1),np.append(Ha,H))
plt.legend(['F','G','H'])

#%%
plt.figure()
lam=np.linalg.eig(J);
plt.plot(np.real(lam[0])*h,np.imag(lam[0])*h,'.');
plt.xlabel('Real(Lambda)')
plt.ylabel('Imag(Lambda)')