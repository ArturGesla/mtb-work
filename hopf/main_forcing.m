clc; clear;
mu=-0.01;
om=1;
b=1;
rss=[];

%
r0=0.01;
th0=0;
K=100;
p=0;
t0=0;
tM=[t0];
yM=[r0;th0];
r00=r0;
%

[t,y]=rk4(@subhopfp,1,1/61,mu,b,om,t0,[r0;th0],p);
t0=t(end);
r0=y(1,end);
th0=y(2,end);
tM=[tM;t0];
yM=[yM,[r0;th0]];

for i=1:16
[t,y]=rk4(@subhopfp,1,1/61,mu,b,om,t0,[r0;th0],p);
% subhopfp(t0,[r0;th0],mu,b,om,p)
% 
% subhopfp(t0,[r0;th0],mu,b,om,0)
r1=y(1,end);
p=0-K*(r1-r0);

r00=r0;
t0=t(end);
r0=y(1,end);
th0=y(2,end);

tM=[tM;t0];
yM=[yM,[r0;th0]];
end
%
plot(tM,yM(1,:),"-x");
xlabel("time"); ylabel("r"); title("Stabilisation of UPO with forcing")