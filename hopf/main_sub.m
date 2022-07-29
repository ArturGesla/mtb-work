
clc; clear;
mu=1;
om=1;
b=1;
rss=[];

%%
mu=-0.30:0.01:0.10;
for i=1:length(mu)
[t,y]=ode45(@(t,y)subhopf(t,y,mu(i),b,om),[0,1],[1; 0]);
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

plot(mu,real(sqrt(1-sqrt(4*mu+1))/sqrt(2)),'kx'); title("Subcritical Hopf bifurcation diagram")