%% l
% lam=1.5;
% lam=(2-sqrt(1.76))/0.4;
% mu=1.9;
% mu=1.99;
mu=2.001;
% mu=2.05;
lam=mu;
x=0;y=0; z=lam;
delta=0.8*lam-0.8*2.8+1;
z=(1-sqrt(delta))/0.4;
r=sqrt(-z*(z-lam));
t=3;
x=r*cos(t);y=r*sin(t);

J=[lam-3+z+0.2*(1-z^2), -1/4, x-0.2*z*2*x;
    1/4,lam-3+z+0.2*(1-z^2),  y-0.2*z*2*y;
    -2*x,-2*y,lam-2*z];
evE=eig(J)
expM=[eig(J)];