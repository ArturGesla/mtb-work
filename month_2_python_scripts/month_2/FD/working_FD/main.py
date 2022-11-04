#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:02:57 2021

@author: gesla
"""
from fd_fun import fd
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp

Re=1e-5;
[x,F,G,H,B,J,J_lam]=fd(Re);

#%
plt.plot(x,F);
plt.plot(x,G);
plt.plot(x,H);
plt.legend(['F','G','H'])

#predictor
lam_dot=1;
u_dot=np.linalg.solve(J,-J_lam)

#%%




n=500; #number of experiments

C=np.zeros(n);
R_array=np.zeros(n);

#bools=np.array([], dtype=bool);

for i in range(0,n):
    Re=1+i*1;
    R_array[i]=Re;
    #ig=res_a.sol(x_plot);
    [x,F,G,H,B,J]=fd(Re,B); 
    C[i]=B[-1];
    print("====== Re=",Re," ======")
    #B=np.append(B,res_a.success);
    
#%%
    
plt.plot(R_array,C)

#%%


#%%