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

def fd(Re=1,ig=0):
    #Re=1000;
    
    n=299;
    l=1;
    h=l/(n+1);
    fields=3;
    
    J=np.zeros((fields*n+1,fields*n+1));
    g=np.zeros(fields*n+1);
    J_lam_out=g;
    
    F=np.zeros(n);
    G=np.zeros(n);
    H=np.zeros(n);
    
    y=np.zeros((fields,n));
    #y[1,:]=np.linspace(0,l,n); # ig for G
    #y[:,:]+=1;
    
    
    B=np.zeros(fields*n+1); #last one is C
    
    for j in range(0,fields):
        for i in range(0,n):
            B[fields*i+j]=y[j,i];
    
    B+=ig;
    
    Fa=0.0;
    Ga=0.0;
    Ha=0.0;
    Fb=0.0;
    Gb=1.0;
    Hb=0.0;
    bc=np.zeros(6);
    bc[0]=Fa;
    bc[1]=Ga;
    bc[2]=Ha;
    bc[3]=Fb;
    bc[4]=Gb;
    bc[5]=Hb;
    
    N=fields*n;
    #last=N;
    
    iter_count=0;
    max_iter=20;
    tol=1;
    
    while(iter_count<max_iter and tol>1e-10):
        [B,tol,J_out,J_lam_out]=fd_iter(J,B,g,bc,N,fields,Re,h);
        iter_count+=1;
        print(iter_count,":",tol)
        if (iter_count==max_iter-1): 
            raise ValueError('Not converged.');
            
    for i in range(0,n):
        F[i]=B[fields*i];
        G[i]=B[fields*i+1];
        H[i]=B[fields*i+2];
        
    F=np.append(np.append(Fa,F),Fb);
    G=np.append(np.append(Ga,G),Gb);
    H=np.append(np.append(Ha,H),Hb);
    x=np.linspace(0,l,n+2);
    
    return [x,F,G,H,B,J_out,J_lam_out];

#%%
    
def fd_iter(J,B,g,bc,N,fields,Re,h):
    last=N;
    C=B[last];
    
    [Fa,Ga,Ha,Fb,Gb,Hb]=bc;
    
    J_lam=np.zeros(len(g));
    
    for i in range(fields,N-fields,fields):        
        g[i]=       (B[i+3]-2*B[i]+B[i-3])/h/h      -Re*(B[i+2]*(B[i+3]-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
        g[i+1]=     (B[i+4]-2*B[i+1]+B[i-2])/h/h    -Re*(B[i+2]*(B[i+4]-B[i-2])/2/h     +2*B[i]*B[i+1]);
        g[i+2]=     (B[i+5]-B[i-1])/2/h             +2*(B[i]+B[i+3])/2;
        
        J_lam[i]=   -(B[i+2]*(B[i+3]-B[i-3])/2/h    +B[i]**2-B[i+1]**2+C);
        J_lam[i+1]= -(B[i+2]*(B[i+4]-B[i-2])/2/h    +2*B[i]*B[i+1]);
        J_lam[i+2]= 0;
        
    
    #Boundaries
    i=0;
    g[i]=       (B[i+3]-2*B[i]+Fa)/h/h      -Re*(B[i+2]*(B[i+3]-Fa)/2/h     +B[i]**2-B[i+1]**2+C);
    g[i+1]=     (B[i+4]-2*B[i+1]+Ga)/h/h    -Re*(B[i+2]*(B[i+4]-Ga)/2/h     +2*B[i]*B[i+1]);
    g[i+2]=     (B[i+5]-Ha)/2/h             +2*(B[i]+B[i+3])/2;
    
    J_lam[i]=   -(B[i+2]*(B[i+3]-Fa)/2/h    +B[i]**2-B[i+1]**2+C);
    J_lam[i+1]= -(B[i+2]*(B[i+4]-Ga)/2/h    +2*B[i]*B[i+1]);
    J_lam[i+2]= 0;

    i=N-fields;
    g[i]=       (Fb-2*B[i]+B[i-3])/h/h      -Re*(B[i+2]*(Fb-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
    g[i+1]=     (Gb-2*B[i+1]+B[i-2])/h/h    -Re*(B[i+2]*(Gb-B[i-2])/2/h     +2*B[i]*B[i+1]);
    g[i+2]=     (Hb-B[i-1])/2/h             +2*(B[i]+Fb)/2;
    
    J_lam[i]=   -(B[i+2]*(Fb-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
    J_lam[i+1]= -(B[i+2]*(Gb-B[i-2])/2/h     +2*B[i]*B[i+1]);
    J_lam[i+2]= 0;
    
    g[last]=    (B[2]-Ha)/h             +2*(B[0]+Fa)/2;
    J_lam[last]=    0;
    
    
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
        J[i+2,i+3]=1;
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
    J[i+2,i+3]=1;
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
    
    J[last,0]=1;
    J[last,2]=1/h;
    #%
    
    dB=-np.linalg.inv(J)@g
    B=B+dB;
    #print(np.linalg.norm(dB))
    return [B,np.linalg.norm(dB),J,J_lam];

def eval_jacobian_u(u,lam,h,bc):
    fields=3;
    N=len(u)-1;
    Re=lam;
    J=np.zeros((len(u),len(u)));
    last=N;
    B=u;
    [Fa,Ga,Ha,Fb,Gb,Hb]=bc;

    
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
        J[i+2,i+3]=1;
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
    J[i+2,i+3]=1;
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
    
    J[last,0]=1;
    J[last,2]=1/h;
    return J;

def eval_jacobian_lam(u,lam,h,bc):
    fields=3;
    N=len(u)-1;
    #Re=lam;
    #J=np.zeros((len(u),len(u)));
    last=N;
    B=u;
    [Fa,Ga,Ha,Fb,Gb,Hb]=bc;

    
    J_lam=np.zeros(len(u));
    
    C=B[last];
    
    for i in range(fields,N-fields,fields):        
#        g[i]=       (B[i+3]-2*B[i]+B[i-3])/h/h      -Re*(B[i+2]*(B[i+3]-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
#        g[i+1]=     (B[i+4]-2*B[i+1]+B[i-2])/h/h    -Re*(B[i+2]*(B[i+4]-B[i-2])/2/h     +2*B[i]*B[i+1]);
#        g[i+2]=     (B[i+5]-B[i-1])/2/h             +2*(B[i]+B[i+3])/2;
        
        J_lam[i]=   -(B[i+2]*(B[i+3]-B[i-3])/2/h    +B[i]**2-B[i+1]**2+C);
        J_lam[i+1]= -(B[i+2]*(B[i+4]-B[i-2])/2/h    +2*B[i]*B[i+1]);
        J_lam[i+2]= 0;
        
    
    #Boundaries
    i=0;
#    g[i]=       (B[i+3]-2*B[i]+Fa)/h/h      -Re*(B[i+2]*(B[i+3]-Fa)/2/h     +B[i]**2-B[i+1]**2+C);
#    g[i+1]=     (B[i+4]-2*B[i+1]+Ga)/h/h    -Re*(B[i+2]*(B[i+4]-Ga)/2/h     +2*B[i]*B[i+1]);
#    g[i+2]=     (B[i+5]-Ha)/2/h             +2*(B[i]+B[i+3])/2;
    
    J_lam[i]=   -(B[i+2]*(B[i+3]-Fa)/2/h    +B[i]**2-B[i+1]**2+C);
    J_lam[i+1]= -(B[i+2]*(B[i+4]-Ga)/2/h    +2*B[i]*B[i+1]);
    J_lam[i+2]= 0;

    i=N-fields;
#    g[i]=       (Fb-2*B[i]+B[i-3])/h/h      -Re*(B[i+2]*(Fb-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
#    g[i+1]=     (Gb-2*B[i+1]+B[i-2])/h/h    -Re*(B[i+2]*(Gb-B[i-2])/2/h     +2*B[i]*B[i+1]);
#    g[i+2]=     (Hb-B[i-1])/2/h             +2*(B[i]+Fb)/2;
    
    J_lam[i]=   -(B[i+2]*(Fb-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
    J_lam[i+1]= -(B[i+2]*(Gb-B[i-2])/2/h     +2*B[i]*B[i+1]);
    J_lam[i+2]= 0;
    
#    g[last]=    (B[2]-Ha)/h             +2*(B[0]+Fa)/2;
    J_lam[last]=    0;
    
    return J_lam;

def eval_G(u,lam,h,bc):
    fields=3;
    N=len(u)-1;
    Re=lam;
    #J=np.zeros((len(u),len(u)));
    last=N;
    B=u;
    [Fa,Ga,Ha,Fb,Gb,Hb]=bc;

    
    #J_lam=np.zeros(len(u));
    g=np.zeros(len(u));
    
    C=B[last];
    
    for i in range(fields,N-fields,fields):        
        g[i]=       (B[i+3]-2*B[i]+B[i-3])/h/h      -Re*(B[i+2]*(B[i+3]-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
        g[i+1]=     (B[i+4]-2*B[i+1]+B[i-2])/h/h    -Re*(B[i+2]*(B[i+4]-B[i-2])/2/h     +2*B[i]*B[i+1]);
        g[i+2]=     (B[i+5]-B[i-1])/2/h             +2*(B[i]+B[i+3])/2;
        
#        J_lam[i]=   -(B[i+2]*(B[i+3]-B[i-3])/2/h    +B[i]**2-B[i+1]**2+C);
#        J_lam[i+1]= -(B[i+2]*(B[i+4]-B[i-2])/2/h    +2*B[i]*B[i+1]);
#        J_lam[i+2]= 0;
        
    
    #Boundaries
    i=0;
    g[i]=       (B[i+3]-2*B[i]+Fa)/h/h      -Re*(B[i+2]*(B[i+3]-Fa)/2/h     +B[i]**2-B[i+1]**2+C);
    g[i+1]=     (B[i+4]-2*B[i+1]+Ga)/h/h    -Re*(B[i+2]*(B[i+4]-Ga)/2/h     +2*B[i]*B[i+1]);
    g[i+2]=     (B[i+5]-Ha)/2/h             +2*(B[i]+B[i+3])/2;
    
#    J_lam[i]=   -(B[i+2]*(B[i+3]-Fa)/2/h    +B[i]**2-B[i+1]**2+C);
#    J_lam[i+1]= -(B[i+2]*(B[i+4]-Ga)/2/h    +2*B[i]*B[i+1]);
#    J_lam[i+2]= 0;

    i=N-fields;
    g[i]=       (Fb-2*B[i]+B[i-3])/h/h      -Re*(B[i+2]*(Fb-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
    g[i+1]=     (Gb-2*B[i+1]+B[i-2])/h/h    -Re*(B[i+2]*(Gb-B[i-2])/2/h     +2*B[i]*B[i+1]);
    g[i+2]=     (Hb-B[i-1])/2/h             +2*(B[i]+Fb)/2;
    
#    J_lam[i]=   -(B[i+2]*(Fb-B[i-3])/2/h     +B[i]**2-B[i+1]**2+C);
#    J_lam[i+1]= -(B[i+2]*(Gb-B[i-2])/2/h     +2*B[i]*B[i+1]);
#    J_lam[i+2]= 0;
    
    g[last]=    (B[2]-Ha)/h             +2*(B[0]+Fa)/2;
#    J_lam[last]=    0;
    
    return g;

def eval_N(u,lam,u_0,lam_0,u_dot,lam_dot,ds):
    N=np.dot((u-u_0),u_dot)+np.dot((lam-lam_0),lam_dot)-ds;
    return N;
    
def assembly(J_u,J_lam,u_dot,lam_dot,G,N):
    n=len(u_dot)+1;
    LHS=np.zeros((n,n));
    RHS=np.zeros(n);
    for i in range(0,n-1):
        for j in range(0,n-1):
            LHS[i,j]=J_u[i,j];
        LHS[i,n-1]=J_lam[i];
        RHS[i]=-G[i];
    RHS[n-1]=-N;
    for j in range(0,n-1):
        LHS[n-1,j]=u_dot[j];
    LHS[n-1,n-1]=lam_dot;
    return [LHS,RHS];

def assembly_vec(J_u,J_lam,u_dot,lam_dot):
    n=len(u_dot)+1;
    LHS=np.zeros((n,n));
    RHS=np.zeros(n);
    for i in range(0,n-1):
        for j in range(0,n-1):
            LHS[i,j]=J_u[i,j];
        LHS[i,n-1]=J_lam[i];
        #RHS[i]=-G[i];
    RHS[n-1]=1;
    for j in range(0,n-1):
        LHS[n-1,j]=u_dot[j];
    LHS[n-1,n-1]=lam_dot;
    return [LHS,RHS];
    
def pc_iteration(u_0,u_dot,lam_0,lam_dot,h,bc,ds):
    #%% predictor
    #ds=-1;
    #ig for u at new lambda
    u=u_0+ds*u_dot; 
    lam=lam_0+ds*lam_dot;
    print("====== Re=",lam,"======")
    
    #%% corrector
    res=1; it=0;
    while (res>1e-6 and it<15):
        J_u=eval_jacobian_u(u,lam,h,bc);
        J_lam=eval_jacobian_lam(u,lam,h,bc);
        G=eval_G(u,lam,h,bc);
        N=eval_N(u,lam,u_0,lam_0,u_dot,lam_dot,ds);
        
        #master matrix
        [LHS,RHS]=assembly(J_u,J_lam,u_dot,lam_dot,G,N);
        dx=np.linalg.solve(LHS,RHS);
        du=dx[0:-1];
        dlam=dx[-1];
        u+=du;
        lam+=dlam;
        res=np.linalg.norm(dx); print(it, ":",res);
        it+=1;
        
    #print("====== Re_fin=",lam,"======")
        
    #% gradient 
    [LHS,RHS]=assembly_vec(J_u,J_lam,u_dot,lam_dot);
    dx=np.linalg.solve(LHS,RHS);
    u_dot=dx[0:-1];
    lam_dot=dx[-1];
    
    #increment
    u_0=u;
    lam_0=lam;
    print("lam_dot= ",lam_dot)
    return [u_0,u_dot,lam_0,lam_dot];

def plot_prof(u,bc,fields,l):
    [Fa,Ga,Ha,Fb,Gb,Hb]=bc;
    
    n=len(u);
    
    F=np.zeros((n-1)//fields);
    G=np.zeros((n-1)//fields);
    H=np.zeros((n-1)//fields);
    
    for i in range(0,(n-1)//fields):
        F[i]=u[fields*i];
        G[i]=u[fields*i+1];
        H[i]=u[fields*i+2];
        
    F=np.append(np.append(Fa,F),Fb);
    G=np.append(np.append(Ga,G),Gb);
    H=np.append(np.append(Ha,H),Hb);
    x=np.linspace(0,l,(n-1)//fields+2);
    
    plt.plot(x,F);
    plt.plot(x,G);
    plt.plot(x,H);
    plt.legend(['F','G','H'])
    
#%%
