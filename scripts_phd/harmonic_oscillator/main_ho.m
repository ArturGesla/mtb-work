clc; clear;
neq=3;
mu=0.04;
om=1;
b=0.5; % aka gamma
rss=[];

m=1;
k=1;
c=1;

%
% r0=0.3;%0.335710687019729+1e-3;
r0=0.0100508962005208-1e-2; %0.142887535036224

th0=0;
K=0.08;
t0=0;
tM=[t0];
yM=[r0;th0];
r00=r0;

T=2*pi/(1-b*mu);
om=sqrt(k/m);
T=2*pi/om;
% dt=5.646783412549308/100;
dt=T/100;
beta=pi/4;

x0=r0;
y0=0;
z0=0;
%
vG=zeros(neq,1)+[0;1/m;0];
[t,y]=rk4(@ho,5*100,dt,m,c,k,t0,[x0;y0;z0],vG);
t0=t(end);
x0=y(1,end);
y0=y(2,end);
tM=[t];
yM=[y];
GM=yM*0;


%
close all;
figure("Position",[2300,600,800,800])
hold on;
plot(tM,yM(1,:),"-");
plot(tM,yM(2,:),"-");
% plot(tM,yM(3,:).^0.5,"-");
% plot(tM,yM(3,:)-mu,"-");

%% pseudo
A=[0,1;-k/m,-c/m];
B=eye(2);

sigma=1i*[0:0.01:5];
res=[];
for i=1:length(sigma)
res=[res;cond(A-sigma(i)*B)];
end

close all;
plot(imag(sigma),res)

%%
close all;
plot3(yM(1,:),yM(2,:),yM(3,:));
grid on;

%%
plot(yM(1,:),yM(2,:));
grid on;
axis square;
axis equal;