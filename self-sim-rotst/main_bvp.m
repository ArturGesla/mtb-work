% xmesh = linspace(0,1,50);
close all;
xmesh = zc;
solinit = bvpinit(xmesh, @guess);
sol = bvp4c(@bvpfcn, @bcfcn, solinit);

plot(sol.x, sol.y(1:3,:), '-o')

f=interp1(sol.x(1,:),sol.y(1,:),zc);
g=interp1(sol.x(1,:),sol.y(2,:),zc);
h=interp1(sol.x(1,:),sol.y(3,:),zw,[],'extrap');
k=mean(sol.y(6,:));
u=zeros(length(f)*4+1,1);
p=gradient(gradient(h,mean(diff(zw))),mean(diff(zw)))/1000-h.*gradient(h,mean(diff(zw)));
u(1:4:end-1)=p;
u(2:4:end)=f;
u(3:4:end)=g;
u(4:4:end)=h;
u(end)=k;



function dydx = bvpfcn(x,y) % equation to solve
dydx = zeros(6,1);
f=y(1); g=y(2); h=y(3);
fp=y(4); gp=y(5); 
%hp=y(6);
%p=y(7); 
k=y(6);
Re=500;
dydx = [fp
    gp
    -2*f
    Re*(k+f^2-g^2+h*fp)
    Re*(2*f*g+h*gp)
    0];
end
%--------------------------------
function res = bcfcn(ya,yb) % boundary conditions
res = [ ya(1)
        ya(2)
        ya(3)
        yb(1)
        yb(2)-1
        yb(3)];
end
%--------------------------------
function g = guess(x) % initial guess for y and y'
g = [0
     -0.1
     0
     0
     0
     0];
end
%--------------------------------