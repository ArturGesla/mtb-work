clc; clear; close all;
%%
neq=2; mu=-0.05; b=-1; omega=1;

% sub
% mu=-0.1;
p = [-1 0 1 0 mu 0];
r = roots(p);
r2=r(3);
r1=r(5);
r=r1;
om=1;
% b=1;
omEff=om+b*r^2;
f=omEff/2/pi;
T=1/f;
TT=T;
T=1;

%
np=1;
dt=T/np;
g=zeros(2*np+1,1);
gstab=zeros(2*np,1);
J=zeros(2*np+1);
Jstab=zeros(2*np);
% J2=zeros(2*np);
u=zeros(2*np+1,1)+T;
% r=sqrt(-mu)
% r=0.05

phi0=pi/8;
x=r*cos(linspace(phi0+0,phi0+2*pi,np+1)); %x=x+rand(1,length(x)).*x*0.1;
y=r*sin(linspace(phi0+0,phi0+2*pi,np+1)); %y=y+rand(1,length(y)).*y*0.1;
u(1:2:end-1)=x(1:end-1); u(2:2:end-1)=y(1:end-1);
% plot(x,y); axis equal;

uM=[];
uM=[uM,u];
%
% calc J and g
 for i=1:1
    for ip=1:np
        
    x=u(ip*2-1); y=u(ip*2); r=sqrt(x^2+y^2);
    xNext=u(mod(ip*2-1+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);
    xPrev=u(mod(ip*2-1-neq-1,neq*np)+1); yPrev=u(mod(ip*2-neq-1,neq*np)+1);
    T=u(2*np+1); dt=T/np; ds=1/np;
%     T=T; dt=T/np; ds=1/np;
    g(ip*2-1)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2))+(xPrev-xNext)/ds/2;
    g(ip*2)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2))+(yPrev-yNext)/ds/2;
    
    gstab(ip*2-1)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2));
    gstab(ip*2)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2));
    
    J(ip*2-1,ip*2-1)=T*((mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x);
    J(ip*2-1,ip*2)=T*(x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y));
    J(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/ds/2;
    J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=1/ds/2;
  
    J(ip*2,ip*2-1)=T*(y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x));
    J(ip*2,ip*2)=T*((mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y);
    J(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/ds/2;
    J(ip*2,mod(ip*2-neq-1,neq*np)+1)=1/ds/2;

    %last column
    J(ip*2-1,end)=((mu+r^2-r^4)*x-y*(omega+b*r^2));
    J(ip*2,end)=((mu+r^2-r^4)*y+x*(omega+b*r^2));
  
    %jstab
    Jstab(ip*2-1,ip*2-1)=T*((mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x);
    Jstab(ip*2-1,ip*2)=T*(x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y));
% %     J(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/ds/2;
%     J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=1/ds/2;
  
    Jstab(ip*2,ip*2-1)=T*(y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x));
    Jstab(ip*2,ip*2)=T*((mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y);
%     J(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/ds/2;
%     J(ip*2,mod(ip*2-neq-1,neq*np)+1)=1/ds/2;
    end
 end
 ustab=u(1:end-1);
 
 
 %%
 ustab'*(Jstab*ustab-gstab*(gstab'*(Jstab*ustab))/norm(gstab)^2)/norm(ustab)
 eig(1/2*(Jstab+transpose(Jstab)))
 Jstab*gstab
 ustab'*Jstab*gstab
 [v,ev]=eig(Jstab); ev
 
 %%
 n=gstab/norm(gstab);
 P=eye(2)-n*n';
 [v,ev]=eig(P*Jstab)