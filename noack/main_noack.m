clc; clear;
neq=3;
mu=0.04;
% mu=15;
% mu=10;
evsAnal=[0;(-1-sqrt(1-8*mu))/2;(-1+sqrt(1-8*mu))/2;]
om=1;
b=1; % aka gamma
% b=0.01; % aka gamma
rss=[];

%
% r0=0.3;%0.335710687019729+1e-3;
r0=0.0100508962005208-1e-2; %0.142887535036224
r0=sqrt(mu)*1.5;

th0=0;
K=0.08;
t0=0;
tM=[t0];
yM=[r0;th0];
r00=r0;

% T=2*pi/(1-b*mu);
T=2*pi;
% dt=5.646783412549308/100;
dt=T/101;
% dt=0.0628;
beta=pi/4;
nt=101;


% x0=r0;
% y0=0;
% z0=mu;

v0=v(1:3)
x0=v0(1)+X(1,1)
y0=v0(2)+X(1,2)
z0=v0(3)+X(1,3)

%
vG=zeros(neq,1);
[t,y]=rk4(@noack,nt,dt,mu,b,om,t0,[x0;y0;z0],vG);
t0=t(end);
x0=y(1,end);
y0=y(2,end);
tM=[t];
yM=[y];
GM=yM*0;


%
close all;
hold on;
plot(tM,yM(1,:),"-");
plot(tM,yM(2,:),"-");
plot(tM,yM(3,:).^0.5,"-");
% plot(tM,yM(3,:)-mu,"-");

%
close all;
plot3(yM(1,:),yM(2,:),yM(3,:));
grid on;
% axis equal;
%%
plot(yM(1,:),yM(2,:));
grid on;
axis square;
axis equal;