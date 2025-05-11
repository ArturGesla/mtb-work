clc; clear; close all;
mu=-0.015;
om=1;
b=-10; % aka gamma
rss=[];

%
% r0=0.3;%0.335710687019729+1e-3;

rex=real(sqrt((-1+sqrt(1+4*mu))/-2));

r0=0.1; %rex*0.1; %0.142887535036224
th0=0;
K=0.04;%sqrt(0.5^2+0.1^2);%0.04;
t0=0;
tM=[t0];
yM=[r0;th0];
r00=r0;

% T=2*pi/(1-b*mu);
T=2*pi/(1+b*sqrt((-1+sqrt(1+4*mu))/-2)^2);
% dt=5.646783412549308/100;
np=100;
dt=T/np;
beta=pi/2;

x0=r0;
y0=0;
%
vG=zeros(2,1);
[t,y]=rk4(@subhopfp,np*1,dt,mu,b,om,t0,[x0;y0],vG);
t0=t(end);
x0=y(1,end);
y0=y(2,end);
tM=[t];
yM=[y];
GM=yM*0;
delay=np;
display(num2str(delay*dt/T))
% plot(sqrt(yM(1,:).^2+yM(2,:).^2))

% period 1 
G=-K*[cos(beta), -sin(beta); sin(beta), cos(beta)];
tic
for i=1:np*50 %10 periods
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
%%
mu=-0.02;
rex2=real(sqrt((-1+sqrt(1+4*mu))/-2));


T2=2*pi/(1+b*sqrt((-1+sqrt(1+4*mu))/-2)^2);
% delay=round(T2/dt);

delay1=delay;
delay2=delay;

dt1=dt;
dt2=1.1*dt;

%% period 2 

tic
for i=1:np*20 %10 periods
[t,y]=rk4(@subhopfp,1,dt1,mu,b,om,t0,[x0;y0],vG);

t0=t(end);
x0=y(1,end);
y0=y(2,end);

tM=[tM,t0];
yM=[yM,[x0;y0]];


%forcing
vG=G*[yM(1,end)-yM(1,end-delay1); yM(2,end)-yM(2,end-delay1)];

GM=[GM,vG];
end
toc

f1=norm(GM(:,end));

tic
for i=1:np*20 %10 periods
[t,y]=rk4(@subhopfp,1,dt2,mu,b,om,t0,[x0;y0],vG);

t0=t(end);
x0=y(1,end);
y0=y(2,end);

tM=[tM,t0];
yM=[yM,[x0;y0]];


%forcing
vG=G*[yM(1,end)-yM(1,end-delay2); yM(2,end)-yM(2,end-delay2)];

GM=[GM,vG];
end
toc

f2=norm(GM(:,end));

%%

%%new delay
%  double a = (alpha1.real() - alpha0.real()) / (lam1 - lam0);
%                     lam = lam0 - alpha0.real() / a;
%                     lam0 = lam1;
%                     alpha0 = alpha1;
%                     lam1 = lam;

% %%delay
% alpha=(f2-f1)/(delay2-delay1);
% delay=delay1-f1/alpha;
% delay1=delay2;
% f1=f2;
% delay2=round(delay);

%%dt
alpha=(f2-f1)/(dt2-dt1);
dt=dt1-f1/alpha;
dt1=dt2;
f1=f2;
dt2=dt;


tic
for i=1:np*20 %10 periods
[t,y]=rk4(@subhopfp,1,dt2,mu,b,om,t0,[x0;y0],vG);

t0=t(end);
x0=y(1,end);
y0=y(2,end);

tM=[tM,t0];
yM=[yM,[x0;y0]];


%forcing
vG=G*[yM(1,end)-yM(1,end-delay2); yM(2,end)-yM(2,end-delay2)];

GM=[GM,vG];
end
toc

f2=norm(GM(:,end));



% minimise the resid wrt T
% hold on; grid on;
% 
% plot(tM/(delay*dt),vecnorm(yM-yM(:,end),1),"-");
% plot(tM(end-np-round(np/5):end-np+round(np/5))/(delay*dt),vecnorm(yM(:,end-np-round(np/5):end-np+round(np/5))-yM(:,end),1),"-");
%
% s=vecnorm(yM(:,end-delay-round(delay/5):end-delay+round(delay/5))-yM(:,end),1);
% [a,c]=min(s);
% delay=round((tM(end)-tM(end-np-round(np/5)+b-1))/dt);
% delay=delay+1;
%%
close all;
hold on; grid on;

plot(tM,vecnorm(yM-yM(:,end),1),"-");
plot(tM(end-delay-round(delay/5):end-delay+round(delay/5)),vecnorm(yM(:,end-delay-round(delay/5):end-delay+round(delay/5))-yM(:,end),1),"-");

s=vecnorm(yM(:,end-delay-round(delay/5):end-delay+round(delay/5))-yM(:,end),1);
[a,c]=min(s);
tM(end)-tM(end-delay-round(delay/5)+c-1)

%% main plot
close all;
f=figure("Position",[2000 200 800 600]);
hold on;
plot(tM/T,yM(1,:),"-");
plot(tM/T,yM(2,:),"-");
plot(tM/T,(tM-tM+1)*rex,"LineWidth",2)
plot(tM/T,(tM-tM+1)*rex2,"LineWidth",2)
xlabel("time in periods"); ylabel("r"); title("UPO stab | K="+num2str(K)+" tau = "+num2str(delay*dt/pi/2,"%4.2f")+" periods | mu="+num2str(mu)+" x0 = "+num2str(r0));  grid on;
r=sqrt(yM(1,:).^2+yM(2,:).^2);
plot(tM/T,r)
legend("x","y","r_{mu=-0.015}","r_{mu=-0.02}");
r(end)-rex
exportgraphics(gcf,"plot.png",'Resolution',150)

%% lsa
close all;
semilogy(tM/T,r); hold on;
% semilogy(tM,rex-r); hold on; 
semilogy(tM/T,r(end)*exp(mu*(t-t(end))))

%% phase
plot(yM(1,:),yM(2,:)); axis equal;
%% force
close all;
% plot(tM/T,GM(1,:),"-");
% plot(tM/T,GM(2,:),"-");
semilogy(tM/T,vecnorm(GM,2),"-"); hold on;
% legend("x force","y force");
title("feedback force"); grid on;
xlabel("time in periods"); ylabel("||force||"); 
exportgraphics(gcf,"plot.png",'Resolution',150)

%%
plot3(yM(1,:),yM(2,:),tM)
%%
r=sqrt(yM(1,:).^2+yM(2,:).^2);
plot(r-rex)
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