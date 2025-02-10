clf;  close all; clear;
uarr=[]; dxarr=[];
%
np=50%*2*2*2*2*2*2;
%%
% xe=0:1/np:1; xe=xe.^1;
% xc=[(-xe(2)+3*xe(1))/2,(xe(1:end-1)+xe(2:end))/2,(3*xe(end)-xe(end-1))/2];
% ! better uni
% dx=1/np; xc=-dx/2:dx:1+dx/2; %harder refinement!
% dx=1/np; xc=0:dx:1; xc=(-cos(xc*pi)+1)/2;k=1; %no nonuni!
dx=1/np; xc=0:dx:1; k=1;
f=@(x) cos(1*pi*x*2*k);
A=zeros(length(xc));
b=zeros(length(xc),1);

%
%  scnd_order;

%fourth_order;
%
% u=A\b;


%pade
fourth_pade; 
u=A\b; u=u(1:2:end); 

%
hold on;
plot(xc,u,'-o'); grid on;

plot(xc,-1/(4*pi*pi)*(f(xc)-1),'-x');
%%
uarr=[uarr;u(np/2+1)];
dxarr=[dxarr;dx];
np=np*2
%%
or=log((uarr(1:end-2)-uarr(2:end-1))./((uarr(2:end-1)-uarr(3:end))))./log(mean(dxarr(1:end-1)./dxarr(2:end)))
%%
clf;
semilogy(abs(or-4),'-x')
%%
plot(uarr)