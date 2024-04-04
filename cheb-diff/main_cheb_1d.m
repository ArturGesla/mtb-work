clc; clear;
uarr=[];
%%
nt=21;
nx=nt;
% x=linspace(-1,1,nx)';
% x=x.^3;
x=-cos(linspace(0,pi,nx))';
%d2udx2-1=0
%
g=x*0;
u=x*0;
J=zeros(length(x),length(x));
% f=x*0-1*x.^2+1;
% f=x*0-1;
f=cos(pi*x);

%
for i=1+1:length(x)-1
    xc=x(i);
    for ik=0:nt-1
%         g(i)=g(i)+u(ik+1)*dTn2(xc,ik)-f(i)*(ik==0);
              g(i)=g(i)+u(ik+1)*dTn2(xc,ik);
        J(i,ik+1)=J(i,ik+1)+dTn2(xc,ik);
    end
    g(i)=g(i)-f(i);
end

i=1;
xc=x(i);
for ik=0:nt-1
    g(i)=g(i)+u(ik+1)*Tn(xc,ik);
    J(i,ik+1)=J(i,ik+1)+Tn(xc,ik);
end

i=length(x);
xc=x(i);
for ik=0:nt-1
    g(i)=g(i)+u(ik+1)*Tn(xc,ik);
    J(i,ik+1)=J(i,ik+1)+Tn(xc,ik);
end
%

u=u-J\g;
% uarr=[uarr;u((length(x)+1)/2)]
uPhys=u*0;
for ik=0:nt-1
    uPhys=uPhys+u(ik+1)*Tn(x,ik);
end
uarr=[uarr;uPhys((length(x)+1)/2)]
%%
d=diff(uarr)
d(1)/d(2)
%%
close all;
plot(x,uPhys,'-x')