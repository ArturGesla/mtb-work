x=[1,2,3,0,0];
y=[4,5,6,0,0];
conv(x,y)
%%
% fftshift
(ifft(fft(x).*fft(y)))
%works but gl to gen it
%%
clc; clear;
rarr=[];
% narr=2.^[1:20];
narr=2.^[1:18];
%
% for n=1:8000;
% for n=narr;
for n=1
% x=rand(1,n);
% y=rand(1,n);

x=[0:4];
y=[5:9];

tic;
res1=conv(x,y);
t1=toc;

% tic;

x=[x,zeros(1,n-1)];
y=[y,zeros(1,n-1)];
tic;
res2=(ifft(fft(x).*fft(y)));
% res2=fft(x);
t2=toc;

rarr=[rarr;[n,t1,t2]];

end
%%
loglog(rarr(:,1),rarr(:,2:end))
% loglog(rarr(:,1),rarr(:,3),'x-');
hold on;
loglog(rarr(:,1),2e-8*narr.*log(narr),'-')