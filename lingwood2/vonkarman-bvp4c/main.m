%%
xmesh = linspace(0,20,100);
solinit = bvpinit(xmesh, @guess);

sol = bvp4c(@bvpfcn, @bcfcn, solinit);

% plot(sol.x, sol.y, '-')
plot(sol.y(1:3,:), sol.x, '-')


function dydx = bvpfcn(x,y) % equation to solve
dydx = zeros(5,1);
dydx = [y(4)
       y(5)
       -2*y(1)
       y(1)^2-y(2)^2+y(4)*y(3)
       2*y(1)*y(2)+y(5)*y(3)];
end
%--------------------------------
function res = bcfcn(ya,yb) % boundary conditions
res = [ya(1)
    ya(2)-1
    ya(3)
    yb(1)
    yb(2) ];
end
%--------------------------------
function g = guess(x) % initial guess for y and y'
g = 1*[x*0
     x*0
     x*0-1
     x*0
     x*0];
end
%-------------------------