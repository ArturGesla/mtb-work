
clc; clear; close all;
lam=1.2;
x=(lam-1)/lam
%
[t,y]=ode45(@(t,y)fun(t,y,lam),[0:0.1:10000],[x+1e-13]);
%
plot(t,y,'-x')
y(end)
y(end)-y(end-1)
%%
z=fft(y-mean(y));
plot(abs((z)))
% xlim([0 1000])
