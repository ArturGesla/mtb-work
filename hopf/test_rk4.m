
clc; clear;
mu=1;
om=1;
b=1;
rss=[];

%%
% mu=-0.30:0.01:0.10;
for i=1:length(mu)
[t,y]=ode45(@(t,y)subhopf(t,y,mu(i),b,om),[0,1],[1; 0]);
r=y(:,1);
rss=[rss;r(end)];
end

for i=1:length(mu)
[tn,yn]=rk4(@subhopf,61,1/61,mu(i),b,om,0,[1;0]);
r=y(:,1);
rss=[rss;r(end)];
end

%%
hold on
plot(t,y(:,1),'-o')
plot(tn,yn(1,:),'x-')