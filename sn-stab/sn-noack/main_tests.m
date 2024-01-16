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
An=@(n) 2*(sigma./n.^2/pi^2/4)./(1+sigma^2./n.^2/pi^2/4).*(exp(sigma)*(1).^n-1);
Bn=@(n) An(n)*2.*n*pi/(-sigma);



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
%%
semilogy(An(1:40),'-x'); hold on; semilogy(-Bn(1:40),'-x'); semilogy(real(z(2:end))*2,'-x'); semilogy(imag(z(2:end))*2,'-sq'); xlim([0,40])

%%

x=0:0.01:2-0.01; y=sin(2*x*pi); ya=exp(x);
z=fftshift(fft(y))/length(y);
za=fftshift(fft(ya))/length(ya);
% z2=conv(z,z,'same');
z2=conv(z,za,'same');
y2=ifft(fftshift(z2))*length(z2); norm(imag(y2))
y2=real(y2);

%%
plot(x,[y;y.^2;y2; y.*ya]'); grid on;


