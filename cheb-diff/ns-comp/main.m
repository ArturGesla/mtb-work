% system("rotst -mode 3 -nx 302 -nz 82 -lx 1 -lz 1 -contipath 1 20 ");
% system("rotst -mode 3 -nx 150 -nz 40 -contipath 9 20 100 200 500 900 1000 1500 2000 2500");
% system("rotst -mode 3 -nx 75 -nz 20 -contipath 4 20 100 200 500");
% system("rotst -mode 3 -nx 76 -nz 21 -contipath 4 20 100 200 500 -uniCosMixedGrid 2 1");
system("rotst -mode 3 -nx 153 -nz 43 -contipath 9 20 100 200 500 900 1000 1500 2000 2500 -uniCosMixedGrid 2 1");
% system("rotst -mode 3 -nx 600 -nz 160 -contipath 4 20 100 200 500");
% system("rotst -mode 100 -re 20 -nx 16 -lx 1 -lz 1  ");

%%
clc; clear;
% a=importdata("u20-600-160.dat");  a0=loadrs2(a,600,160);
% a=importdata("u500-600-160.dat");  a0=loadrs2(a,600,160);
% a=importdata("u20-302-82.dat");  a0=loadrs2(a,302,82);
% a=importdata("u500-300-80.dat");  a1=loadrs2(a,300,80);
a=importdata("u2500-150-40.dat");  a2=loadrs2(a,150,40);
a=importdata("u2500-153-43.dat");  a2m=loadrs2(a,153,43);
% a=importdata("u500-75-20.dat");  a3=loadrs2(a,75,20);
% a=importdata("u500-76-21.dat");  a3m=loadrs2(a,76,21);

%%
b=importdata("u-10.dat"); n=10; n=n;
bx=reshape(b(:,1),[n+1,n]); bx=(bx+1)/2;
bz=reshape(b(:,2),[n+1,n]); bz=(bz+1)/2;
bu=reshape(b(:,3),[n+1,n]);

%%
clf;
b2=importdata("u-12-uni.dat"); n=12; n=n;
b2x=reshape(b2(:,1),[n,n]*10); b2x=(b2x+1)/2;
b2z=reshape(b2(:,2),[n,n]*10); b2z=(b2z+1)/2;
b2u=reshape(b2(:,3),[n,n]*10);
contour(b2x,b2z,b2u,[-1:0.2:1]/20,'r'); hold on;


b2=importdata("u-10-uni.dat"); n=10; n=n;
b2x=reshape(b2(:,1),[n,n]*10); b2x=(b2x+1)/2;
b2z=reshape(b2(:,2),[n,n]*10); b2z=(b2z+1)/2;
b2u=reshape(b2(:,3),[n,n]*10);
% contour(b2x,b2z,b2u,[-1:0.2:1]/20,'k');


b2=importdata("u-6-uni.dat"); n=6; n=n;
b2x=reshape(b2(:,1),[n,n]*10); b2x=(b2x+1)/2;
b2z=reshape(b2(:,2),[n,n]*10); b2z=(b2z+1)/2;
b2u=reshape(b2(:,3),[n,n]*10);
% contour(b2x,b2z,b2u,[-1:0.2:1]/20,'g');


b2=importdata("u-16-uni.dat"); n=16; n=n;
b2x=reshape(b2(:,1),[n,n]*10); b2x=(b2x+1)/2;
b2z=reshape(b2(:,2),[n,n]*10); b2z=(b2z+1)/2;
b2u=reshape(b2(:,3),[n,n]*10);
contour(b2x,b2z,b2u,[-1:0.2:1]/20,'c');
%%
clf; 
% a=a3; mesh(a.xc,a.zc,a.u); hold on;
a=a3m; mesh(a.xc,a.zc,a.u); hold on;
%  mesh(bx,bz,bu);
%  mesh(b2x,b2z,b2u);
%%
clf; hold on; 
% grid on;
% a=a0; contour(a.xc,a.zc,a.v,[0:1:10],'c')
% a=a1; contour(a.xc,a.zc,a.v,[0:1:10],'r')
% a=a2; contour(a.xc,a.zc,a.v,[0:1:10],'k')
% a=a3; contour(a.xc,a.zc,a.v,[0:1:10],'b')
% a=a0; contour(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1),[-1:0.2:1],'c')
% a=a1; contour(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1),[-1:0.2:1],'r')
% a=a2; contour(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1),[-1:0.2:1],'k')
% a=a0; contour(a.xu(1,:),a.zc(:,1),a.u(:,1:end-1),[-1:0.2:1]/20,'c')
% a=a0; contour(bx,bz,bu,[-1:0.2:1]/20,'r')
% a=a0; contour(b2x,b2z,b2u,[-1:0.2:1]/20,'r')

% a=a0; contour(a.xu(1,:),a.zc(:,1),a.u(:,1:end-1),[-1:0.2:1],'r')
% a=a1; contour(a.xu(1,:),a.zc(:,1),a.u(:,1:end-1),[-1:0.2:1],'r')
% a=a2; contour(a.xu(1,:),a.zc(:,1),a.u(:,1:end-1),[-1:0.2:1],'g')
% a=a3; contour(a.xu(1,:),a.zc(:,1),a.u(:,1:end-1),[-1:0.2:1],'b')
% a=a2m; contour(a.xu(1,:),a.zc(:,1),a.u(:,1:end-1),[-1:0.2:1],'k')

% a=a0; contour(a.xc(1,:),a.zw(:,1),a.w(1:end-1,:),[-0.2:0.02:0.1],'r')
% a=a1; contour(a.xc(1,:),a.zw(:,1),a.w(1:end-1,:),[-0.2:0.02:0.1],'r')
a=a2; contour(a.xc(1,:),a.zw(:,1),a.w(1:end-1,:),[-0.2:0.04:0.1],'k')

mesh(a.xc,a.zc,a.zc-100,'EdgeColor','k')
% a=a3; contour(a.xc(1,:),a.zw(:,1),a.w(1:end-1,:),[-0.2:0.02:0.1],'b')
% a=a3m; contour(a.xc(1,:),a.zw(:,1),a.w(1:end-1,:),[-0.2:0.02:0.1],'k')
a=a2m; contour(a.xc(1,:),a.zw(:,1),a.w(1:end-1,:),[-0.2:0.04:0.1],'g')
% mesh(a.xc,a.zc,a.zc-100,'EdgeColor','k')
% a=a2; contour(a.xu,a.zc,a.u(2:end-1,2:end-1),[0:1:10],'k')
%%
a=a3;
pcolor(a.xc,a.zc,a.u);  colormap(parula(8))
%%
plot(a.w(:,1:end/10:end),'x-')
