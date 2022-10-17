
clc; clear;
mu=1;
om=1;
b=1;
rss=[];

%%
mu=-0.30:0.01:0.10;
for i=1:length(mu)
[t,y]=ode45(@(t,y)subhopf(t,y,mu(i),b,om),[0,100],[1; 0]);
r=y(:,1);
rss=[rss;r(end)];
end
rssleft=rss; rss=[];
for i=1:length(mu)
[t,y]=ode45(@(t,y)subhopf(t,y,mu(i),b,om),[0,1000],[1e-2; 0]);
r=y(:,1);
rss=[rss;r(end)];
end
rssright=rss;
%%
plot(mu,rss,'x-')
%%
hold on; grid on; grid minor;
plot(mu,rssleft,'x'); xlabel("mu"); ylabel("r")
plot(mu,rssright,'x')
rssu=(sqrt(1-sqrt(4*mu+1))/sqrt(2));

plot(mu,real(sqrt(1-sqrt(4*mu+1))/sqrt(2)),'kx'); title("Subcritical Hopf bifurcation diagram")
%%
close all;

f=figure(Position=[2200 202 911 598]); fnts=14;
hold on;
set(f,'defaulttextinterpreter','latex')


hold on; grid on; grid minor;
plot(mu(6:end),rssleft(6:end),'b-x',"LineWidth",2 ); xlabel("$\mu$"); ylabel("r")
plot(mu(1:end-10),rssright(1:end-10),'-bx',"LineWidth",2 );

plot(mu(6:end),real(rssu(6:end)),'rx--',"LineWidth",2 ); title("Subcritical Hopf bifurcation diagram $\vert$ $\omega=1$ $\vert$ $b=1$")

set(gca,"FontSize",fnts,"FontName","Latin Modern Math");

exportgraphics(gcf,'plot.png','Resolution',300)