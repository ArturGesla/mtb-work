function dydt=odefun(t,y,mu,b,om)
dydt=zeros(2,1);
dydt=[mu*y(1)-y(1)^3; om+b*y(1)^2];
end