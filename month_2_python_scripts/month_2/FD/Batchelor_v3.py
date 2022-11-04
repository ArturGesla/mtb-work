# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 12:25:01 2021

@author: Artur
"""
#%%
import numpy as np
import matplotlib.pyplot as plt

Re=100;

n=99;
l=1;
h=l/(n+1);
fields=6;

J=np.zeros((fields*n,fields*n));
g=np.zeros(fields*n);

#F=np.zeros(n);
#G=np.zeros(n)+1;
#H=np.zeros(n);

y=np.zeros((fields,n));
y[1,:]+=1; # ig for G

B=np.zeros(fields*n);

for j in range(0,fields):
    for i in range(0,n):
        B[fields*i+j]=y[j,i];

Fa=0.0;
Ga=0.0;
Ha=0.0;
Fb=0.0;
Gb=1.0;
Hb=0.0;

#%%
N=fields*n;
for i in range(fields,N-fields,fields):
    #for j in range(0,fields):
    #    g[i+j]=(B[i+j+fields]-B[i+j-fields])/2/h;
    
    #g[i]+=      -B[i+3];
    #g[i+1]+=    -B[i+4];
    #g[i+2]+=    2*B[i];
    #g[i+3]+=    -Re*(B[i]**2-B[i+1]**2+B[i+5]+B[i+2]*B[i+3]);
    #g[i+4]+=    -Re*(2*B[i]*B[i+1]+B[i+4]*B[i+2]);
    #g[i+5]+=    0;
        
    g[i]=       (B[i+6]-B[i-6])/2/h         -B[i+3];
    g[i+1]=     (B[i+7]-B[i-5])/2/h         -B[i+4];
    g[i+2]=     (B[i+8]-B[i-4])/2/h         +2*B[i];
    g[i+3]=     (B[i+9]-B[i-3])/2/h         -Re*(B[i]**2-B[i+1]**2+B[i+5]+B[i+2]*B[i+3]);
    g[i+4]=     (B[i+10]-B[i-2])/2/h        -Re*(2*B[i]*B[i+1]+B[i+4]*B[i+2]);
    g[i+5]=     (B[i+11]-B[i-1])/2/h        +0;
        
    

#Boundaries
i=0;
g[i]=       (B[i+6]-Fa)/2/h     -B[i+3];
g[i+1]=     (B[i+7]-Ga)/2/h     -B[i+4];
g[i+2]=     (B[i+8]-Ha)/2/h     +2*B[i];
g[i+3]=     (B[i+9]-B[i+3])/h   -Re*(B[i]**2-B[i+1]**2+B[i+5]+B[i+2]*B[i+3]);
g[i+4]=     (B[i+10]-B[i+4])/h  -Re*(2*B[i]*B[i+1]+B[i+4]*B[i+2]);
g[i+5]=     (B[i+11]-B[i+5])/h  +0;


i=N-fields;
g[i]=       (Fb-B[i-6])/2/h     -B[i+3];
g[i+1]=     (Gb-B[i-5])/2/h     -B[i+4];
g[i+2]=     (Hb-B[i-4])/2/h     +2*B[i];
g[i+3]=     (B[i+3]-B[i-3])/h   -Re*(B[i]**2-B[i+1]**2+B[i+5]+B[i+2]*B[i+3]);
g[i+4]=     (B[i+4]-B[i-2])/h   -Re*(2*B[i]*B[i+1]+B[i+4]*B[i+2]);
g[i+5]=     (B[i+5]-B[i-1])/h   +0;



#%
for i in range(fields,N-fields,fields):
    J[i,i+6]=1/2/h;         J[i,i-6]=-1/2/h;        J[i,i+3]=-1;
    J[i+1,i+7]=1/2/h;       J[i+1,i-5]=-1/2/h;      J[i+1,i+4]=-1;
    J[i+2,i+8]=1/2/h;       J[i+2,i-4]=-1/2/h;      J[i+2,i]=2;
    J[i+3,i+9]=1/2/h;       J[i+3,i-3]=-1/2/h;      J[i+3,i]=-Re*(2*B[i]);      J[i+3,i+1]=Re*(2*B[i+1]);   J[i+3,i+2]=-Re*(B[i+3]);  J[i+3,i+3]=-Re*(B[i+2]);  J[i+3,i+5]=-Re;
    J[i+4,i+10]=1/2/h;      J[i+4,i-2]=-1/2/h;      J[i+4,i]=-Re*(2*B[i+1]);    J[i+4,i+1]=-Re*(2*B[i]);    J[i+4,i+2]=-Re*(B[i+4]);  J[i+4,i+4]=-Re*(B[i+2]);  
    J[i+5,i+11]=1/2/h;      J[i+5,i-1]=-1/2/h;      

#Boundaries
i=0;
J[i,i+6]=1/2/h;                                 J[i,i+3]=-1;
J[i+1,i+7]=1/2/h;                               J[i+1,i+4]=-1;
J[i+2,i+8]=1/2/h;                               J[i+2,i]=2;
J[i+3,i+9]=1/h;       J[i+3,i+3]=-1/h;          J[i+3,i]=-Re*(2*B[i]);      J[i+3,i+1]=Re*(2*B[i+1]);   J[i+3,i+2]=-Re*(B[i+3]);  J[i+3,i+3]=-Re*(B[i+2]);  J[i+3,i+5]=-Re;
J[i+4,i+10]=1/h;      J[i+4,i+4]=-1/h;          J[i+4,i]=-Re*(2*B[i+1]);    J[i+4,i+1]=-Re*(2*B[i]);    J[i+4,i+2]=-Re*(B[i+4]);  J[i+4,i+4]=-Re*(B[i+2]);  
J[i+5,i+11]=1/h;      J[i+5,i+5]=-1/h;      

i=N-fields;
J[i,i-6]=-1/2/h;        J[i,i+3]=-1;
J[i+1,i-5]=-1/2/h;      J[i+1,i+4]=-1;
J[i+2,i-4]=-1/2/h;      J[i+2,i]=2;

J[i+3,i+3]=1/h;         J[i+3,i-3]=-1/h;        J[i+3,i]=-Re*(2*B[i]);      J[i+3,i+1]=Re*(2*B[i+1]);   J[i+3,i+2]=-Re*(B[i+3]);  J[i+3,i+3]=-Re*(B[i+2]);  J[i+3,i+5]=-Re;
J[i+4,i+4]=1/h;         J[i+4,i-2]=-1/h;        J[i+4,i]=-Re*(2*B[i+1]);    J[i+4,i+1]=-Re*(2*B[i]);    J[i+4,i+2]=-Re*(B[i+4]);  J[i+4,i+4]=-Re*(B[i+2]);  
J[i+5,i+5]=1/h;         J[i+5,i-1]=-1/h;      
    

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