x=0:1/10:1-1/10; y=exp(0*x).*(sin(1*x*2*pi)+0.2*sin(2*x*2*pi)+0.1*sin(3*x*2*pi));
z=fft(y)/length(y);
% z=fft(y);
close all; plot([0:length(x)-1],abs(z(1:end)),'x-'); 
% hold on; plot(imag(z),'x-');


y=exp(1*x).*(sin(1*x*2*pi)+0.2*sin(2*x*2*pi)+0.1*sin(3*x*2*pi));
z=fft(y)/length(y);
% z=fft(y);
hold on; plot([0:length(x)-1],abs(z(1:end)),'x-'); 

%%
close all;
plot(real(z(2:end)))

z2=[z(1:end/2+1),zeros(1,40),z(end/2+1:end)]
x2=length(z2)*ifft(z2);

plot([0:1/(length(x2)-1):1],x2,'x-');
hold on;
plot(x,y,'-o')

%%
