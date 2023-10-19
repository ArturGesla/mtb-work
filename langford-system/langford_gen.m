
a=-3; b=-9.3; c=8; d=-3; e=5.98; close all;
langfordG = @(t,y) [a*y(1)+b*y(2)+y(1)*y(3);
    c*y(1)+d*y(2)+y(2)*y(3);
    e*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))]; % Anonymous Function
opt=odeset('RelTol',1e-10);
 [t,X] = ode45(langford, [0,100], [ -4,0,3],opt);
% plot(X(:,1),X(:,2))
plot3(X(:,1),X(:,2),X(:,3))

%%
% close all;
nt=30000; %generalised does not work wtf xd
dt=1e-1;
t0=0;
% a=-3; b=-9.3; c=8; d=-3; e=5.98; close all;
% a=-3; b=-1; c=-b; d=a; e=2/3; close all; x0=[1e-2;0;1e-3]; %per orbit
 mu=2.01; a=mu-3; b=-1/1; c=-b; d=a; e=mu; c2=0.2; close all; x0=[0.1;0;2]; %per orbit
% a=-3; b=-8; c=8; d=-3; e=5.98; close all; x0=[0.1;0.1;0.1]; %per orbit
langfordG = @(t,y) [a*y(1)+b*y(2)+y(1)*y(3);
    c*y(1)+d*y(2)+y(2)*y(3);
    e*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))];
langfordG = @(t,y) [a*y(1)+b*y(2)+y(1)*y(3)+c2*y(1)*(1-y(3)^2);
                    c*y(1)+d*y(2)+y(2)*y(3)+c2*y(2)*(1-y(3)^2);
                    e*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))];% Anonymous Function
% [t,y]=rk4_2(@(t,y)langfordG(t,y),nt,dt,t0,[1;1;1]);
[t,y]=rk4_2(@(t,y)langfordG(t,y),nt,dt,t0,x0);
% [t,y]=rk4_2(@(t,y)langfordG(t,y),nt,dt,t0,[0.003; 0.001; 0.001]);
% plot3(y(1,:),y(2,:),y(3,:));
% xlabel("x"); ylabel("y");  zlabel("z"); 
%
subplot(1,2,1)
plot(t,y'); hold on;
plot(t,sqrt(y(1,:).^2+y(2,:).^2));
subplot(1,2,2)
plot3(y(1,:),y(2,:),y(3,:));
delta=0.8*mu-0.8*2.8+1
z=[-1-sqrt(delta),-1+sqrt(delta)]/-0.4
r=sqrt(-z.*(z-mu))
sqrt(y(1,end).^2+y(2,end).^2)
y(3,end)
%%
plot3(y(1,:),y(2,:),y(3,:));
%%
close all;
figure("Position",[1921           1        1920         997])
for i=1:10
subplot(2,5,i)
% l=2.01; 
% l=1.98+i/10*(2.023-1.98); 
l=1.5+i/10*(1.98-1.5); 
langford = @(t,y) [(l-3)*y(1)-1/4*y(2)+y(1)*(y(3)+0.2*(1-y(3)^2));
    (l-3)*y(2)+1/4*y(1)+y(2)*(y(3)+0.2*(1-y(3)^2));
    l*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))]; % Anonymous Function
 [t,X] = ode45(langford, [0,1000], [ 0.2334 -0.9401 1.0858]);
% plot(X(:,1),X(:,2))
plot3(X(:,1),X(:,2),X(:,3))
xlim([-1.2 1.2]); ylim([-1.2 1.2]); zlim([0 2]);
title("lam: "+num2str(l));
grid on;
end
%%
exportgraphics(gcf,"plot1.png","Resolution",300);

%%
X=y';
x=X(end-20000:end,1);
y=X(end-20000:end,2);
z1=X(end-20000:end,3);
t1=t(end-20000:end);
% plot3(x,y,z)
plot(t1,x); clf;
[f,z]=ft(t1,x); semilogy(f,z); hold on; grid on;
[f,z]=ft(t1,y); semilogy(f,z); hold on; grid on;
[f,z]=ft(t1,z1); semilogy(f,z); hold on; grid on;
% xlim([0 0.2]);
% ylim([1e-6 1])


%% stab
lam=1.5;
lam=(2-sqrt(1.76))/0.4;
lam=2.01;
x=0;y=0; z=lam;
delta=0.8*lam-0.8*2.8+1;
z=(1-sqrt(delta))/0.4;
r=sqrt(-z*(z-lam));
t=4;
x=r*cos(t);y=r*sin(t);

J=[lam-3+z+0.2*(1-z^2), -1/4, x-0.2*z*2*x;
    1/4,lam-3+z+0.2*(1-z^2),  y-0.2*z*2*y;
    -2*x,-2*y,lam-2*z];
eig(J)