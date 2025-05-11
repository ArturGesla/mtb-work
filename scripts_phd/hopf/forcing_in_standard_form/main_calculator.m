%% sub
mu=-0.015;
p = [-1 0 1 0 mu 0];
r = roots(p);
r=r(end);
om=1;
b=-10;
omEff=om+b*r^2;
f=omEff/2/pi;
T=1/f;