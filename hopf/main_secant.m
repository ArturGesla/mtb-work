%%
muc=-0.1;
r0=0;
r1=1;

[t0,y0]=ode45(@(t,y)subhopf(t,y,muc,b,om),[0,1000],[r0; 0]);
[t1,y1]=ode45(@(t,y)subhopf(t,y,muc,b,om),[0,1000],[r1; 0]);
r0=y0(end,1);
r1=y1(end,1);
alpha=0.5;
rh=[r0;r1;];

for i=1:100
r=alpha*r0+(1-alpha)*r1; rh=[rh;r];
[t,y]=ode45(@(t,y)subhopf(t,y,muc,b,om),[0,1000],[r; 0]);
rk=y(end,1);
if(abs(rk-r1)>abs(rk-r0))
    r0=r; 
else
    r1=r;
end
end
%%
hold on;
plot(t0,y0(:,1),'-x')
plot(t1,y1(:,1),'-x')
%%
plot(rh); xlabel("iteration of secant method"); ylabel("r"); title("History of r"); grid on;
%%
semilogy(abs(diff(rh))); xlabel("iteration of secant method"); ylabel("dr"); title("History of r"); grid on;
%%
plot(t,y(:,1)); xlabel("time"); ylabel("r"); grid on; title("r history ")
%%
plot(t,y(:,1),"-x"); xlabel("time"); ylabel("r"); grid on; title("r history ")