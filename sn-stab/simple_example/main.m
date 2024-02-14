
nt=20;
t=0:2*pi/nt:2*pi; t=t';
xex=1+cos(t)*0.9;%1+2*cos(t)+1/2+1/2*cos(2*t); 
% xex=xex.^2;

plot(t,xex);

% xex=xex*0;
% fft
z=fft(xex(1:end-1))./(length(xex)-1);

n=5;
u=[real(z(1:n));imag(z(1:n))];
%%

% [g,jac]=calculateRhsAndJacLinSN(1,n,u,1);
[g,jac]=calculateRhsAndJacNLSN(1,n,u,1);
u=u-jac\g;
norm(g)
%%

