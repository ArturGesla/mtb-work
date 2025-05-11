clc; clear;
mu=-0.02;
om=1;
b=-10; % aka gamma
rss=[];

%
% r0=0.3;%0.335710687019729+1e-3;
r0=0.100508962005208-1e-2; %0.142887535036224

th0=0;
K=0.08;
t0=0;
tM=[t0];
yM=[r0;th0];
r00=r0;

T=2*pi/(1-b*mu);
% dt=5.646783412549308/100;
dt=T/100;
beta=pi/4;

x0=r0;
y0=0;
%
vG=zeros(2,1);
[t,y]=rk4(@subhopfp_r3,403-1,dt,mu,b,om,t0,[x0;y0],vG);
t0=t(end);
x0=y(1,end);
y0=y(2,end);
tM=[t];
yM=[y];
GM=yM*0;

% plot(sqrt(yM(1,:).^2+yM(2,:).^2))

G=-K*[cos(beta), -sin(beta); sin(beta), cos(beta)];

for i=1:10000
[t,y]=rk4(@subhopfp_r3,1,dt,mu,b,om,t0,[x0;y0],vG);

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
plot(tM,yM(1,:),"-");
plot(tM,yM(2,:),"-");
xlabel("time"); ylabel("r"); title("UPO stab | K="+num2str(K)+" mu="+num2str(mu)+" x0 = "+num2str(r0)); legend("x","y"); grid on;
%%
hold on;
plot(GM(1,:),"-");
plot(GM(2,:),"-");
legend("x force","y force");
title("feedback force"); grid on;
%%
plot3(yM(1,:),yM(2,:),tM)
%%
r=sqrt(yM(1,:).^2+yM(2,:).^2);
plot(r)
%%
r(end)-sqrt(-mu)
norm(GM(:,end))