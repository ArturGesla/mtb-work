#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 08:39:07 2021

@author: gesla
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

x1=2;
y1=1;
x2=1;
y2=2;
l1=-1;
l2=-10;

a11=(y2*x1*l1-l2*y1*x2)/(x1*y2-y1*x2);
a12=(x1*x2*(l2-l1))/(x1*y2-y1*x2);
a21=(y1*y2*(l1-l2))/(x1*y2-y1*x2);
a22=(y2*x1*l2-l1*y1*x2)/(x1*y2-y1*x2);

A=np.zeros((2,2));
A[0,0]=a11;
A[0,1]=a12;
A[1,0]=a21;
A[1,1]=a22;

lam=np.linalg.eig(A)
l=lam[0];
v=lam[1];
v1=v[:,0];
v2=v[:,1];
l1=l[0];
l2=l[1];

#%%

t=np.arange(0,10,0.01);
x=np.zeros((2,len(t)));
for i in range(0,len(t)):
    x[:,i]=v1*np.exp(t[i]*l1)-v2*np.exp(t[i]*l2);    
    
#%%
plt.plot(t,x[1,:])
plt.plot(t,x[0,:])
    
#%%
plt.plot(x[0,:],x[1,:])
    
    