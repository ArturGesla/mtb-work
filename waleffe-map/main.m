clc; clear;
mu=0.6; eps=2;
xleps=((1-mu)/2^eps)^(1/(1+eps)); xreps=1-xleps;

%%
n=1000000; %10 mln takes a minute
x0=rand([n,1])*(xreps-xleps)+xleps;
x=x0;
%% integration
% x=x0; %hold on;
%x0=0.2
% s1=~(x0>0.5);
% s2=x0>0.5;
for i=1:40
x0=(~(x0>0.5)).*x0.*mu./(1-2^eps.*(x0).^(1+eps))+(x0>0.5).*(1-x0).*mu./(1-2^eps.*(1-x0).^(1+eps));
% plot(x0)
x=[x,x0];
display(i);
end
%%
lt=(xleps<x)&(x<xreps); lt=sum(lt,2);
histogram(lt,30); set(gca,'yscale','log'); grid on;
title("Histogram of lifetimes"); ylabel("number of samples"); xlabel("number of iterations")

%%
%  xhat=
k=0;
slope=(x(:,3:end)-x(:,2:end-1))./(x(:,2:end-1)-x(:,1:end-2))+k*(x(:,2:end-1)-x(:,1:end-2));
xhat=(x(:,2:end-1)-slope.*x(:,1:end-2))./(1-slope);
%%
histogram(xhat(:,:),linspace(0.1,1,1000)); grid on; grid minor;

%%
syms x;
a=solve(mu*(1-x)/(1-2^eps*(1-x)^(1+eps))-x==0,'Maxdegree',4);
a=eval(a)
b=solve(mu*(x)/(1-2^eps*x^(1+eps))-x==0)
b=eval(b)
x1=a(2); x2=b(2);
%%
histogram(xhat(:,:),linspace(0.4,0.6,121)); grid on; grid minor;
hold on; plot([x1,x1],[0,n/2]); plot([x2,x2],[0,n/2]); legend("hist of xhat","x1="+num2str(x1),"x2="+num2str(x2)); 
xlabel("x"); ylabel("number of samples"); title("Histogram of xhat") 

