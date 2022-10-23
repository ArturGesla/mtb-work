
clc; clear; close all;
lam=3.1;
x=(lam-1)/lam
%
% [t,y]=ode45(@(t,y)fun(t,y,lam),[0:0.1:1000],[x+1e-13]);
y=[x+1e-13]
for  i=1:1000
    y(end+1)=y(end)*lam*(1-y(end));
end
%
plot(y,'-x')
y(end)
y(end)-y(end-1)
