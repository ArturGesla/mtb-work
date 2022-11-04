#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:36:04 2021

@author: gesla
"""
 
from fd_fun import fd
import matplotlib.pyplot as plt
import numpy as np
from scipy.linalg import eigh 

Re=0.01;
[x,F,G,H,B,J]=fd(Re);

#%
plt.plot(x,F);
plt.plot(x,G);
plt.plot(x,H);
plt.legend(['F','G','H'])

RHS=np.zeros(J.shape);
RHS[range(0,len(RHS),3),range(0,len(RHS),3)]=1;
RHS[range(1,len(RHS),3),range(1,len(RHS),3)]=1;
RHS[-1,-1]=0;

lam=eigh(J,RHS);
plt.plot(np.real(lam[0]),np.imag(lam[0]),'x');
plt.xlabel('Real(Lambda)')
plt.ylabel('Imag(Lambda)')