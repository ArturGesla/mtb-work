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
np=100;
dt=T/np;
g=zeros(2*np,1);
J=zeros(2*np);
J2=zeros(2*np);
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
 for i=1:16
    for ip=1:np-1
    x=u(ip*2-1); y=u(ip*2); r=sqrt(x^2+y^2);
    xNext=u(mod(ip*2-1+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);
    xPrev=u(mod(ip*2-1-neq-1,neq*np)+1); yPrev=u(mod(ip*2-neq-1,neq*np)+1);
    xPrevPrev=u(mod(ip*2-1-neq*2-1,neq*np)+1); yPrevPrev=u(mod(ip*2-neq*2-1,neq*np)+1);

    g(ip*2-1)=(mu+r^2-r^4)*x-y*(omega+b*r^2)-(3*x-4*xPrev+xPrevPrev)/dt/2;
    g(ip*2)=(mu+r^2-r^4)*y+x*(omega+b*r^2)-(3*y-4*yPrev+yPrevPrev)/dt/2;
    
    J(ip*2-1,ip*2-1)=(mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x;
    J(ip*2-1,ip*2)=x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y);
    J(ip*2-1,ip*2-1)=-3/dt/2;
    J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=4/dt/2;
    J(ip*2-1,mod(ip*2-1-neq*2-1,neq*np)+1)=-1/dt/2;
  
    J(ip*2,ip*2-1)=y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x);
    J(ip*2,ip*2)=(mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y;
    J(ip*2,ip*2)=-3/dt/2;
    J(ip*2,mod(ip*2-neq-1,neq*np)+1)=4/dt/2;
    J(ip*2,mod(ip*2-neq*2-1,neq*np)+1)=-1/dt/2;

%     J2(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/dt/2;
%     J2(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=1/dt/2;
%   
%     J2(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/dt/2;
%     J2(ip*2,mod(ip*2-neq-1,neq*np)+1)=1/dt/2;
  
    end
    ip=np;
     x=u(ip*2-1); y=u(ip*2); r=sqrt(x^2+y^2);
    xNext=u(mod(ip*2-1+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);
    xPrev=u(mod(ip*2-1-neq-1,neq*np)+1); yPrev=u(mod(ip*2-neq-1,neq*np)+1);
    xPrevPrev=u(mod(ip*2-1-neq*2-1,neq*np)+1); yPrevPrev=u(mod(ip*2-neq*2-1,neq*np)+1);

    g(ip*2-1)=(mu+r^2-r^4)*x-y*(omega+b*r^2)-(3*x-4*xPrev+xPrevPrev)/dt/2;
%     g(ip*2)=(mu+r^2-r^4)*y+x*(omega+b*r^2)-(3*y-4*yPrev+yPrevPrev)/dt/2;
    g(ip*2)=-(3*y-4*yPrev+yPrevPrev)/dt/2;
    
    J(ip*2-1,ip*2-1)=(mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x;
    J(ip*2-1,ip*2)=x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y);
    J(ip*2-1,ip*2-1)=-3/dt/2;
    J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=4/dt/2;
    J(ip*2-1,mod(ip*2-1-neq*2-1,neq*np)+1)=-1/dt/2;
  
%     J(ip*2,ip*2-1)=y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x);
%     J(ip*2,ip*2)=1;
    J(ip*2,ip*2)=-3/dt/2;
    J(ip*2,mod(ip*2-neq-1,neq*np)+1)=4/dt/2;
    J(ip*2,mod(ip*2-neq*2-1,neq*np)+1)=-1/dt/2;
    %
    du=-J\g;
   
    norm(du)
    u=u+du;
    uM=[uM,u];
 end
%%
hold on;
% u=u-du2;
plot(u(1:2:end),u(2:2:end),'x-')
% plot(u(1:2:end)+utang(1:2:end),u(2:2:end)+utang(2:2:end))
% plot([u(1:2:end), u(1:2:end)+utang(1:2:end)]',[u(2:2:end),u(2:2:end)+utang(2:2:end)]'); axis equal
plot([u(1:2:end), u(1:2:end)+du(1:2:end)]',[u(2:2:end),u(2:2:end)+du(2:2:end)]','-o'); axis equal
% plot([u(1:2:end), u(1:2:end)+du(1:2:end)]',[u(2:2:end),u(2:2:end)+du(2:2:end)]','-x'); 

%%
plot(uM(1:2:end,:),uM(2:2:end,:))
%%
close all;
f=figure(Position=[2200 202 911 598]); fnts=14;
% f=figure(Position=[2200 202 300 300]); fnts=14;
% hold on;
set(f,'defaulttextinterpreter','latex')


hold on;
plot(uM(1:2:end,1),uM(2:2:end,1))
plot(uM(1:2:end,end),uM(2:2:end,end),'-o')

x=r*cos(linspace(0,2*pi,np+1));
y=r*sin(linspace(0,2*pi,np+1));
u(1:2:end)=x(1:end-1); u(2:2:end)=y(1:end-1);
plot(u(1:2:end),u(2:2:end),'-x')
axis equal;
% grid on;

h1=legend("iteration=0","iteration=8, resid=1e-16","analytical solution",Location="best");
set(h1, 'Interpreter','latex')
grid on; grid minor; 
set(gca,"FontSize",fnts,"FontName","Latin Modern Math");
xlabel("x"); ylabel("y"); title("Multi timestep newton method $ \vert$ subHoph $\mu=-0.015 $ $ \vert$ 100 points for period")
exportgraphics(gcf,'plot.png','Resolution',200)

%%
hold on;
plot(uM(1,:)',uM(2,:)','-x'); axis equal;
plot(u(1:2:end),u(2:2:end))
%%
plot(sqrt(u(1:2:end).^2+u(2:2:end).^2))
%%
hold on;
plot(uM')
plot(sqrt(sum(uM.^2,1))')
