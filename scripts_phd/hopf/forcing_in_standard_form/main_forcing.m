clc; clear;
mu=-0.1;
om=1;
b=1;
rss=[];

%
r0=0.335710687019729;
th0=0;
K=20;
G=zeros(2);
t0=0;
tM=[t0];
yM=[r0;th0];
r00=r0;
dt=5.646783412549308/100;
%

[t,y]=rk4(@subhopfp,103-1,dt,mu,b,om,t0,[r0;th0],G);
t0=t(end);
r0=y(1,end);
th0=y(2,end);
tM=[t];
yM=[y];
GM=yM*0;

for i=1:10
[t,y]=rk4(@subhopfp,1,dt,mu,b,om,t0,[r0;th0],G);

t0=t(end);
r0=y(1,end);
th0=y(2,end);

tM=[tM,t0];
yM=[yM,[r0;th0]];

%forcing
G(1,1)=-K*(yM(1,end)-yM(1,end-100));
ph1=yM(2,end)-yM(2,end-90);
ph2=ph1-2*pi;
G(2,2)=-K*(ph2);

GM=[GM,[G(1,1);G(2,2)]];
end
%
hold on;
plot(yM(1,:),"-x");
plot(yM(2,:)/2/pi,"-x");
xlabel("time"); ylabel("r"); title("Stabilisation of UPO with forcing")
%%
figure(); 
plot(GM(1,:),"-");
title("r force");
figure(); 
plot(GM(2,:),"-");
title("th force")
