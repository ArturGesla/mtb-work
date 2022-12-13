mu=1;
om=1;
b=1;
rss=[];

%
mu=-0.105:0.01:0.105;
for i=1:length(mu)
[t,y]=ode45(@(t,y)superhopf(t,y,mu(i),b,om),[0,1000],[0.1; 0]);
r=y(:,1);
rss=[rss;r(end)];
end
%%
plot(mu,rss,'x-')