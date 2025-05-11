function [t,y]=rk4_2(odefun,nt,dt,t0,y0)
% dydt=odef(1,1,mu,b,om);
t=[t0];
y=y0;

for i=1:nt
k1=odefun(t0,y0);
k2=odefun(t0+dt/2,y0+k1*dt/2 );
k3=odefun(t0+dt/2,y0+k2*dt/2 );
k4=odefun(t0+dt,y0+k3*dt );
y0=y0+1/6*(k1+2*k2+2*k3+k4)*dt;
t0=t0+dt;
t(end+1)=t0;
y(:,end+1)=y0;
end

end