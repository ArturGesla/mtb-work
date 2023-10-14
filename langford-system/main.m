clc; clear;
lam=2.02;
k=1;
c=1;
neq=3;

x0=1e-6;
y0=0;
z0=1e-5;
t0=0;
dt=0.1;
%
vG=zeros(neq,1);
v=[0.2334
   -0.9401
    1.0858];
nt=10000;
% [t,y]=rk4(@langford,nt,dt,lam,c,k,t0,[x0;y0;z0],vG);
[t,y]=rk4(@langford,nt,dt,lam,c,k,t0,v,vG);
% vG=zeros(neq+1,1);
% [t1,y1]=rk4(@langford2,nt,dt,lam,c,k,t0,[v;v(end)^2],vG);
%
clf;
 plot(t,y'); hold on;
% plot(t1,y1','--')