#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 13:36:30 2021

@author: gesla
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 12:25:01 2021

@author: Artur
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

Re=100;

n=299;
l=1;
h=l/(n+1);
fields=3;

J=np.zeros((fields*n+1,fields*n+1));
g=np.zeros(fields*n+1);

F=np.zeros(n);
G=np.zeros(n);
H=np.zeros(n);

y=np.zeros((fields,n));
y[1,:]=np.linspace(0,l,n); # ig for G
#y[:,:]+=1;

B=np.zeros(fields*n+1); #last one is C

for j in range(0,fields):
    for i in range(0,n):
        B[fields*i+j]=y[j,i];

Fa=0.0;
Ga=0.0;
Ha=0.0;
Fb=0.0;
Gb=1.0;
Hb=0.0;

N=fields*n;
last=N;
#%%
C=B[last];
for i in range(fields,N-fields,fields):        
    g[i]=       (B[i+3]-2*B[i]+B[i-3])/h/h      -Re*(B[i+2]*(B[i+3]-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
    g[i+1]=     (B[i+4]-2*B[i+1]+B[i-2])/h/h    -Re*(B[i+2]*(B[i+4]-B[i-2])/2/h     +2*B[i]*B[i+1]);
    g[i+2]=     (B[i+5]-B[i-1])/2/h             +2*(B[i]);

    

#Boundaries
i=0;
g[i]=       (B[i+3]-2*B[i]+Fa)/h/h      -Re*(B[i+2]*(B[i+3]-Fa)/2/h     +B[i]**2-B[i+1]**2+C);
g[i+1]=     (B[i+4]-2*B[i+1]+Ga)/h/h    -Re*(B[i+2]*(B[i+4]-Ga)/2/h     +2*B[i]*B[i+1]);
g[i+2]=     (B[i+5]-Ha)/2/h             +2*(B[i]);


i=N-fields;
g[i]=       (Fb-2*B[i]+B[i-3])/h/h      -Re*(B[i+2]*(Fb-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
g[i+1]=     (Gb-2*B[i+1]+B[i-2])/h/h    -Re*(B[i+2]*(Gb-B[i-2])/2/h     +2*B[i]*B[i+1]);
g[i+2]=     (Hb-B[i-1])/2/h             +2*(B[i]);

g[last]=    (B[2]-Ha)/h             +2*(Fa)/2;



#%
for i in range(fields,N-fields,fields):
    J[i,i]=     -2/h/h      -Re*2*B[i];
    J[i,i+1]=               -Re*(-2*B[i+1]);
    J[i,i+2]=               -Re*((B[i+3]-B[i-3])/2/h);
    J[i,i+3]=   1/h/h       -Re*B[i+2]/2/h;
    J[i,i-3]=   1/h/h       -Re*(-B[i+2]/2/h);
    J[i,last]=              -Re;
    
    J[i+1,i]=               -Re*2*B[i+1];
    J[i+1,i+1]=-2/h/h       -Re*2*B[i];
    J[i+1,i+2]=             -Re*(B[i+4]-B[i-2])/2/h;
    J[i+1,i+4]=1/h/h        -Re*B[i+2]/2/h;
    J[i+1,i-2]=1/h/h        -Re*(-B[i+2]/2/h);
    
    J[i+2,i]=1;
    #J[i+2,i+3]=1;
    J[i+2,i+5]=1/2/h;
    J[i+2,i-1]=-1/2/h;
#Boundaries
i=0;
J[i,i]=     -2/h/h      -Re*2*B[i];
J[i,i+1]=               -Re*(-2*B[i+1]);
J[i,i+2]=               -Re*((B[i+3]-Fa)/2/h);
J[i,i+3]=   1/h/h       -Re*B[i+2]/2/h;
#J[i,i-3]=   1/h/h       -Re*(-B[i+2]/2/h);
J[i,last]=-Re;

J[i+1,i]=               -Re*2*B[i+1];
J[i+1,i+1]=-2/h/h       -Re*2*B[i];
J[i+1,i+2]=             -Re*(B[i+4]-Ga)/2/h;
J[i+1,i+4]=1/h/h        -Re*B[i+2]/2/h;
#J[i+1,i-2]=1/h/h        -Re(-B[i+2]/2/h);

J[i+2,i]=1;
#J[i+2,i+3]=1;
J[i+2,i+5]=1/2/h;
#J[i+2,i-1]=-1/2/h;

i=N-fields;
J[i,i]=     -2/h/h      -Re*2*B[i];
J[i,i+1]=               -Re*(-2*B[i+1]);
J[i,i+2]=               -Re*((Fb-B[i-3])/2/h);
#J[i,i+3]=   1/h/h       -Re*B[i+2]/2/h;
J[i,i-3]=   1/h/h       -Re*(-B[i+2]/2/h);
J[i,last]=              -Re;

J[i+1,i]=               -Re*2*B[i+1];
J[i+1,i+1]=-2/h/h       -Re*2*B[i];
J[i+1,i+2]=             -Re*(Gb-B[i-2])/2/h;
#J[i+1,i+4]=1/h/h        -Re*B[i+2]/2/h;
J[i+1,i-2]=1/h/h        -Re*(-B[i+2]/2/h);

J[i+2,i]=1;
#J[i+2,i+3]=1;
#J[i+2,i+5]=1/2/h;
J[i+2,i-1]=-1/2/h;

#J[last,0]=1;
J[last,2]=1/h;
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
lam=np.linalg.eig(J);
plt.plot(np.real(lam[0])*h,np.imag(lam[0])*h,'x');
plt.xlabel('Real(Lambda)')
plt.ylabel('Imag(Lambda)')