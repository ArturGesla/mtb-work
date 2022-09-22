k=0:0.01:0.04;
gamma=-10;
mu=-0.015;

hold on;

beta=pi/4;
f=1+k*2*pi/(1-gamma*mu)*(cos(beta)+gamma*sin(beta));
plot(k,f)

beta=pi/8*3;
f=1+k*2*pi/(1-gamma*mu)*(cos(beta)+gamma*sin(beta));
plot(k,f)

beta=pi/2;
f=1+k*2*pi/(1-gamma*mu)*(cos(beta)+gamma*sin(beta));
plot(k,f)

legend("pi/4","pi*3/8","pi/2")
grid on;

mu=-0.03;

hold on;

beta=pi/4;
f=1+k*2*pi/(1-gamma*mu)*(cos(beta)+gamma*sin(beta));
plot(k,f)

beta=pi/8*3;
f=1+k*2*pi/(1-gamma*mu)*(cos(beta)+gamma*sin(beta));
plot(k,f)

beta=pi/2;
f=1+k*2*pi/(1-gamma*mu)*(cos(beta)+gamma*sin(beta));
plot(k,f)

legend("pi/4","pi*3/8","pi/2")
grid on;