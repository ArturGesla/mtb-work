clc; clear;
mu=-0.05;
om=1;
b=-1; % aka gamma
rss=[];

%
% r0=0.3;%0.335710687019729+1e-3;
p = [-1 0 1 0 mu 0];
r = roots(p);
r2=r(3);
r1=r(5);

rex=real(sqrt((-1+sqrt(1+4*mu))/-2));

r0=rex+1e-6;%0.1; %rex*0.1; %0.142887535036224
th0=0;
K=0.0;%sqrt(0.5^2+0.1^2);%0.04;
t0=0;
tM=[t0];
yM=[r0;th0];
r00=r0;

% T=2*pi/(1-b*mu);
T=2*pi/(1+b*sqrt((-1+sqrt(1+4*mu))/-2)^2);
% dt=5.646783412549308/100;
np=240;
dt=T/np;
beta=pi/2;

x0=r0;
y0=0;
%
vG=zeros(2,1);
[t,y]=rk4(@subhopfp,np*50,dt,mu,b,om,t0,[x0;y0],vG);
t0=t(end);
x0=y(1,end);
y0=y(2,end);
tM=[t];
yM=[y];
GM=yM*0;
delay=24;
display(num2str(delay*dt/T))
% plot(sqrt(yM(1,:).^2+yM(2,:).^2))
%
G=-K*[cos(beta), -sin(beta); sin(beta), cos(beta)];
tic
for i=1:np*1 %10 periods
[t,y]=rk4(@subhopfp,1,dt,mu,b,om,t0,[x0;y0],vG);

t0=t(end);
x0=y(1,end);
y0=y(2,end);

tM=[tM,t0];
yM=[yM,[x0;y0]];

%forcing
% G(1,1)=-K*(yM(1,end)-yM(1,end-100));
% G(2,2)=-K*(yM(2,end)-yM(2,end-100));

%forcing
vG=G*[yM(1,end)-yM(1,end-delay); yM(2,end)-yM(2,end-delay)];

GM=[GM,vG];
end
toc
%
close all;
f=figure("Position",[20 20 800 600]);
hold on;
plot(tM/T,yM(1,:),"-");
plot(tM/T,yM(2,:),"-");
plot(tM/T,(tM-tM+1)*rex)
xlabel("time in periods"); ylabel("r"); title("UPO stab | K="+num2str(K)+" tau = "+num2str(delay*dt/pi/2,"%4.2f")+" periods | mu="+num2str(mu)+" x0 = "+num2str(r0));  grid on;
r=sqrt(yM(1,:).^2+yM(2,:).^2);
plot(tM/T,r)
legend("x","y","r_{exact}","r");
r(end)-rex
exportgraphics(gcf,"subhopfti.png",'Resolution',100)

%% lsa
close all;
semilogy(tM/T,r); hold on;
% semilogy(tM,rex-r); hold on; 
semilogy(tM/T,r(end)*exp(mu*(t-t(end))))

%% phase
plot(yM(1,:),yM(2,:)); axis equal;
%%
hold on;
plot(tM/T,GM(1,:),"-");
plot(tM/T,GM(2,:),"-");
legend("x force","y force");
title("feedback force"); grid on;
%%
plot3(yM(1,:),yM(2,:),tM)
%%
r=sqrt(yM(1,:).^2+yM(2,:).^2);
% plot(r-rex)
% semilogy(tM,-(r-rex)); hold on;
% semilogy(tM,(r)-r1); hold on; semilogy(tM,1e-6*exp(tM*0.0944)); 
semilogy(tM,r2-(r)); hold on; semilogy(tM,1e10*exp((tM-100)*-1.6944)); 
% set(gca,"Yscale","lin");
ylim([1e-16 1])
r(end)-rex
grid on;
xlabel("time"), ylabel("r-rExact"); title("subHopf mu=-0.05; b=-1; omega=1");
exportgraphics(gcf,"subHopfstab1.png","resolution",150)
%%
r=sqrt(yM(1,:).^2+yM(2,:).^2);
plot(r)
semilogy(r)
r(end)
plot(tM(1000:1500),log(r(1000:1500)))

%%
Teff=dt*1;

lam=(mu+Teff*K*(cos(beta)*mu+sin(beta)))/(1+2*Teff*K*cos(beta)+Teff*K)
%%
K=-20:0.01:10
lam=(mu+Teff*K*(cos(beta)*mu+sin(beta)))./(1+2*Teff*K*cos(beta)+Teff*K)
plot(K,lam)
%%
Teff=dt*1;
A=[mu, -1; 1, mu]; 
B=[1+K*cos(beta)*Teff, -K*sin(beta)*Teff; K*sin(beta)*Teff, 1+K*cos(beta)*Teff];
eig(A,B)