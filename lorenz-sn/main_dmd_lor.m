% dt=0.1;
% t=0:dt:30*pi;
% T=2*pi;
% om=2*pi/T;
% % y=exp(0.156*t);
% ev=1e-2+0.1i;
% y=sin(om*t).*(exp(ev*t));
% z=cos(om*t).*(exp(ev*t));
% 
% z=exp(1i*om*t).*exp(ev*t);
% y=real(z);
% z=imag(z);
% % close all; plot(t,y);
% close all; plot(y,z);
%%
clc; close all; clear;
%%
np=100; r=24;
main_lorenz_ti
%%
y=X';
%%
V1=y(:,1:end-1);
V2=y(:,2:end);
[U,S,V]=svd(V1);
M=U'*V2*V*pinv(S)
eig(M)
% ev=log(eig(M))/dt