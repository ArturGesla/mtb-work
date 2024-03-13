clc; clear;
uarr=[];
%%
nx=101;
x=linspace(-1,1,nx)';
x=x.^3;
% x=-cos(linspace(0,pi,nx))';
%d2udx2-1=0
g=x*0;
u=x*0;
J=zeros(length(x),length(x));
f=x*0-1*x.^2+1;
% f=x*0-1*x+1;
for i=1+1:length(x)-1
    dx1=x(i)-x(i-1);
    dx2=x(i+1)-x(i);
g(i)=-2/dx1/dx2*u(i)+2/dx2/(dx1+dx2)*u(i+1)+2/dx1/(dx1+dx2)*u(i-1)-f(i);
J(i,i)=-2/dx1/dx2;
J(i,i-1)=+2/dx1/(dx1+dx2);
J(i,i+1)=+2/dx2/(dx1+dx2);
end
i=1; g(i)=u(i); J(i,i)=1;
i=length(x); g(i)=u(i); J(i,i)=1;

u=u-J\g;
uarr=[uarr;u((length(x)+1)/2)]
%%
d=diff(uarr)
d(1)/d(2)
%%
close all;
plot(x,u,'-x')