clc; close all;  clear;

%%

dx=0.001;
t=0:dx:1-dx;
sigma=1;



% x=exp(sigma*abs(t));
% A0=1/sigma*(exp(sigma)-1);
% An=@(n) 2*(sigma./n.^2/pi^2)./(1+sigma^2./n.^2/pi^2).*(exp(sigma)*(-1).^n-1);
% Bn=@(n)0;Z

x=exp(sigma*mod(t,1));
A0=1/sigma*(exp(sigma)-1);
An=@(n) 2*(sigma/n^2/pi^2/4)/(1+sigma^2/n^2/pi^2/4)*(exp(sigma)*(1)^n-1);
Bn=@(n) An(n)*2*n*pi/(-sigma);



N=2;
xn=A0;
for i=1:N
xn=xn+An(i)*cos(2*pi*i*t)+Bn(i)*sin(2*pi*i*t);
% xn=xn+An(i)*cos(pi*i*t)+Bn(i)*sin(pi*i*t);
end


clf; 
plot(t,x,'-sq'); 
hold on; grid on;
plot(t,xn,'-+'); 

%

z=fft(x)/length(x);
zc=z;
n=N+1;
zc=z*0; zc([1:n,end-n+2:end])=z([1:n,end-n+2:end]);
xc=ifft(zc)*length(zc);

plot(t,xc,'-x');

% clf; plot(xn-xc)