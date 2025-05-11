% addpath('C:\Users\Artur\Documents\GitHub\rotst2\scripts\source_for_mtb')
addpath('/people/gesla/Documents/git/rotst2/scripts/source_for_mtb')
%%
clf; clc; clear;
n=2; m=10; 
% n=1; m=1;
% n=2; m=4;
% n=10; m=10;
nx=600; nz=160; 
% nx=2*m+1; nz=2*n+1;
ifig=0;

%%
% a00=rand(1); 
% ase=rand(n,m)+rand(n,m)*1i;
% as=rand(n,1)+rand(n,1)*1i;
% ae=rand(1,m)+rand(1,m)*1i;
% ane=rand(n,m)+rand(n,m)*1i;
% rng(ifig)
a00=0; %rand(1); 
ase=exp(1i*rand(n,m)*2*pi);
as=exp(1i*rand(n,1)*2*pi); %rand(n,1)+rand(n,1)*1i;
ae=exp(1i*rand(1,m)*2*pi); % rand(1,m)+rand(1,m)*1i;
ane=exp(1i*rand(n,m)*2*pi); % rand(n,m)+rand(n,m)*1i;

z=[a00,ae,zeros(1,nx-m-m-1),fliplr(conj(ae));
    as,ase,zeros(n,nx-m-m-1),conj(fliplr((ane)));
    zeros(nz-n-n-1,nx);
conj(flipud(as)),(flipud(ane)),zeros(n,nx-m-m-1),conj(fliplr(flipud(ase)))];
x=ifft2(z);
x=x/max(max(x));

xu=linspace(0,10,nx); zw=linspace(0,1,nz);

% %simpler
% [X,Z]=meshgrid(xu,zw);  xs=linspace(0,10,21); zs=linspace(0,1,5); xx=rand(length(zs),length(xs))*2-1; x=interp2(xs,zs,xx,X,Z);
% 
% %cheb
% xcheb=x*0;
% for i=0:5-1
%     for j=0:21-1
%         xcheb=xcheb+xx(i+1,j+1)*cos(j*acos(X/10*2-1)).*cos(i*acos(Z/1*2-1));
%     end
% end
% x=xcheb; x=x/max(max(x));


h=pcolor(xu,zw,x); h.EdgeAlpha=0; %shading interp; 
colormap(parula(8)); c=colorbar(); axis equal; xlim([0 10]); ylim([0,1]);
c.TickLabelInterpreter='latex'; c.Limits=([-1 1]); xlabel('r'); ylabel('z'); 

%
fnts=10; jfm_plt_aid_comm; size_recth15; ifig=ifig+1;

%%
exportgraphics(gcf,"rand_f_field"+num2str(ifig)+".png")

