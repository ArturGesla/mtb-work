function dydt=odefun(t,y,r,b,sigma)
dydt=zeros(3,1);
% lorenz = @(t,x) [10*(x(2)-x(1));r*x(1)-1*x(1)*x(3)-x(2);x(1)*x(2)-(8/3)*x(3)];    % Anonymous Function
dydt=[sigma*(y(2)-y(1));r*y(1)-y(1)*y(3)-y(2);y(1)*y(2)-(b)*y(3)];
end