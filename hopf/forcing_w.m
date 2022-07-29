clc; clear;
mu=-0.1;
om=1;
b=1;
rss=[];

%
r0=0.3;
th0=0;
%controler
w0=0.001;
k=20;
lamc=0.05;
t0=0;
tM=[t0];
yM=[r0;th0;w0];
%
r00=r0;
th00=r0;
w00=w0;

for i=1:4000
[t,y]=rk4w(@subhopfw,1,0.1,mu,b,om,t0,[r0;th0;w0],lamc,k,[r00;th00;w00]);

r00=r0;
th00=r0;
w00=w0;

t0=t(end);
r0=y(1,end);
th0=y(2,end);
w0=y(3,end);

tM=[tM;t0];
yM=[yM,[r0;th0;w0]];
end
%
hold on;
% plot(tM,yM(1,:),"-");
plot(tM,(tM-tM+1)*sqrt(1 - sqrt(1 + 4 *mu))/sqrt(2),"-"); grid on; grid minor;
plot(tM,yM(1,:),"-");
legend("exact","forced","Location","best")
xlabel("time"); ylabel("r"); title("Stabilisation of UPO with forcing nonlinear | mu="+num2str(mu))