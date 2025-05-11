clc; clear; close all;
neq=2; mu=-0.015; b=-10;omega=1;

% sub
% mu=-0.1;
p = [-1 0 1 0 mu 0];
r = roots(p);
r=r(end);
om=1;
% b=1;
omEff=om+b*r^2;
f=omEff/2/pi;
T=1/f;
%
np=300;
dt=T/np;
g=zeros(2*np,1);
J=zeros(2*np);
u=zeros(2*np,1)+3;
% r=sqrt(-mu)
% r=0.05

x=r*cos(linspace(0,2*pi,np+1)); x=x+rand(1,length(x)).*x*0.1;
y=r*sin(linspace(0,2*pi,np+1)); y=y+rand(1,length(y)).*y*0.1;
u(1:2:end)=x(1:end-1); u(2:2:end)=y(1:end-1);
% plot(x,y); axis equal;

uM=[];
uM=[uM,u];
%%
% calc J and g
 for i=1:10
    for ip=1:np
    x=u(ip*2-1); y=u(ip*2); r=sqrt(x^2+y^2);
    xNext=u(mod(ip*2-1+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);
    g(ip*2-1)=(mu+r^2-r^4)*x-y*(omega+b*r^2)+(x-xNext)/dt;
    g(ip*2)=(mu+r^2-r^4)*y+x*(omega+b*r^2)+(y-yNext)/dt;
    
    J(ip*2-1,ip*2-1)=(mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x+1/dt;
    J(ip*2-1,ip*2)=x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y);
    J(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/dt;
  
    J(ip*2,ip*2-1)=y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x);
    J(ip*2,ip*2)=(mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y+1/dt;
    J(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/dt;
    end
    %
    du=-J\g;
    norm(du)
    u=u+du;
    uM=[uM,u];
 end
%%
hold on;
plot(u(1:2:end),u(2:2:end))
%%
plot(uM(1:2:end,:),uM(2:2:end,:))
%%
plot(uM(1,:)',uM(2,:)'); axis equal;
%%
plot(sqrt(u(1:2:end).^2+u(2:2:end).^2))
%%
hold on;
plot(uM')
plot(sqrt(sum(uM.^2,1))')
