%% sub
mu=-0.1;
p = [-1 0 1 0 mu 0];
r = roots(p);
r=r(end);
om=1;
b=1;
omEff=om+b*r^2;
f=omEff/2/pi;