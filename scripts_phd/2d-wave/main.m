
kx=-1; ky=1;
lx=2*pi/abs(kx); ly=2*pi/abs(ky); n=10;
x=linspace(0,lx,n); x=x(1:end-1);
y=linspace(0,ly,n+10); y=y(1:end-1);

[X,Y]=meshgrid(x,y);
Z=sin(X*kx*2+Y*ky);
%%
surf(X,Y,Z)

%%
plot(real(fft(Z(:,:)')))
%%
z=fft2(Z);
%%
surf(abs(z))

%% 1d naive
clc; clear;
cd /home/gesla/Documents/rotst-work/plq;
addpath /home/gesla/Documents/git/rotst2/scripts/source_for_mtb

% a=importdata("u960-256-128.dat"); 
a=importdata("evcReal-256-128-960-0.480000i.dat"); 

a=loadrs2(a,256,128);

vp=a.v;
pcolor(a.xc,a.zc,vp); shading interp; colormap(parula(6)); colorbar; 
%caxis([-max(max(abs(vp(2:end-1,2:end-1))))*0 max(max(abs(vp(2:end-1,2:end-1))))])

delta=0.4
r2=5/(1-delta)
% caxis([0 r2])

%%

tiledlayout(2,1); 
nexttile;
zz=fft(vp(:,:));
zz=zz./vecnorm(zz);
pcolor(log10(abs(zz(1:6,:)))); shading interp; colormap([parula(16)]); colorbar()

nexttile;
% a=loadrs2(a,256,128);

vp=a.v;
pcolor(a.xc,1-a.zc,vp); shading interp; colormap(parula(16)); colorbar; caxis([-1 1]*1e-2)
