clc; clear;
mu=-0.02;
om=1;
b=-10; % aka gamma
rex=real(sqrt((-1+sqrt(1+4*mu))/-2));
T=2*pi/(1+b*sqrt((-1+sqrt(1+4*mu))/-2)^2);

lags=[T T];
k=0.0;
beta=pi/4;

tspan = [0 5*T];
tic
sol = dde23(@(t,y,Z)ddefun(t,y,Z,mu,b,om,k,beta), lags, @history, tspan);
toc
t=sol.x;
x=sol.y(1,:);
y=sol.y(2,:);

r0=0;

hold on;
plot(t/T,x,"-");
plot(t/T,y,"-");
xlabel("time in periods"); ylabel("r"); title("UPO stab | K="+num2str(k)+" mu="+num2str(mu)+" x0 = "+num2str(r0)); legend("x","y"); grid on;






%%


rss=[];

%
% r0=0.3;%0.335710687019729+1e-3;



r0=0.1*rex; %0.142887535036224
th0=0;
K=0.09;%sqrt(0.5^2+0.1^2);%0.04;
t0=0;
tM=[t0];
yM=[r0;th0];
r00=r0;

% T=2*pi/(1-b*mu);
T=2*pi/(1+b*sqrt((-1+sqrt(1+4*mu))/-2)^2);
% dt=5.646783412549308/100;
dt=T/100;
beta=pi/4;

x0=r0;
y0=0;
%
vG=zeros(2,1);
[t,y]=rk4(@subhopfp,103-1,dt,mu,b,om,t0,[x0;y0],vG);
t0=t(end);
x0=y(1,end);
y0=y(2,end);
tM=[t];
yM=[y];
GM=yM*0;

% plot(sqrt(yM(1,:).^2+yM(2,:).^2))

G=-K*[cos(beta), -sin(beta); sin(beta), cos(beta)];

for i=1:30000
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
vG=G*[yM(1,end)-yM(1,end-100); yM(2,end)-yM(2,end-100)];

GM=[GM,vG];
end
%
hold on;
plot(tM/T,yM(1,:),"-");
plot(tM/T,yM(2,:),"-");
xlabel("time in periods"); ylabel("r"); title("UPO stab | K="+num2str(K)+" mu="+num2str(mu)+" x0 = "+num2str(r0)); legend("x","y"); grid on;
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
plot(r)
r(end)-rex
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