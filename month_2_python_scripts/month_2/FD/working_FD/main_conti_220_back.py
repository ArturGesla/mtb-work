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
from fd_fun import eval_jacobian_u,eval_jacobian_lam,eval_G,eval_N,assembly,assembly_vec,pc_iteration,plot_prof
from numpy.linalg import norm

n=99;
l=1;
h=l/(n+1);
fields=3;
size=fields*n+1;

G=np.zeros(size);
G_u=np.zeros((size,size));
G_lambda=np.zeros(size);

#u=np.zeros(size);
u_dot=np.zeros(size);
u_0=np.zeros(size);
#lam=2;
lam_dot=1;
lam_0=500;
N=0;
ds=0;

bc=np.array((0,0,0,0,1,0)); #[Fa,Ga,Ha,Fb,Gb,Hb]

#%% Solve for u_0, u_dot, lam_dot - whole 0th state is known
res=1; it=0;
while (res>1e-6 and it<15):
    J_u=eval_jacobian_u(u_0,lam_0,h,bc);
    J_lam=eval_jacobian_lam(u_0,lam_0,h,bc);
    G=eval_G(u_0,lam_0,h,bc);
    #N=eval_N(u,lam,u_0,lam_0,u_dot,lam_dot,ds);
    
    du=np.linalg.solve(J_u,-G);
    u_0+=du;
    res=norm(du); print(it, ":",res);
    it+=1;

#% gradient 
[LHS,RHS]=assembly_vec(J_u,J_lam,u_dot,lam_dot);
dx=np.linalg.solve(LHS,RHS);
u_dot=dx[0:-1];
lam_dot=dx[-1];

plot_prof(u_0,bc,fields,l);
#%%
for i in range(500,220,-1):
    [u_0,u_dot,lam_0,lam_dot]=pc_iteration(u_0,u_dot,lam_0,lam_dot,h,bc);
    
#plot_prof(u_0,bc,fields,l);

#%%
u_220=          u_0;
u_dot_220=      u_dot;
lam_220=        lam_0;
lam_dot_220=    lam_dot;