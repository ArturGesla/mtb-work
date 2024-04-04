clc; clear;
uarr=[];
%%
nx=81; ny=nx;
% x=linspace(-1,1,nx)'; y=x;
% x=x.^3;
x=-cos(linspace(0,pi,nx))';y=x;

%d2udx2-1=0
g=zeros(nx*ny,1);
u=g;
J=sparse(nx*ny,nx*ny,0);
%
% f=x*0-1*x.^2+1;
% f=x*0-1*x+1;
f=u*0-1;
% f=cos(pi*x);
%
for ix=1+1:length(x)-1
    for iy=1+1:length(y)-1
        dx1=x(ix)-x(ix-1); dx2=x(ix+1)-x(ix); dy1=y(iy)-y(iy-1); dy2=y(iy+1)-x(iy);
ip=iy+(ix-1)*ny;
ipxp=iy+(ix-1+1)*ny;
ipxm=iy+(ix-1-1)*ny;
ipyp=iy+(ix-1)*ny+1;
ipym=iy+(ix-1)*ny-1;

g(ip)=g(ip)+-2/dx1/dx2*u(ip)+2/dx2/(dx1+dx2)*u(ipxp)+2/dx1/(dx1+dx2)*u(ipxm);
g(ip)=g(ip)+-2/dy1/dy2*u(ip)+2/dy2/(dy1+dy2)*u(ipyp)+2/dy1/(dy1+dy2)*u(ipym);
g(ip)=g(ip)+-f(ip);

J(ip,ip)=J(ip,ip)+-2/dx1/dx2-2/dy1/dy2;
J(ip,ipxp)=J(ip,ipxp)+2/dx2/(dx1+dx2);
J(ip,ipxm)=J(ip,ipxm)+2/dx1/(dx1+dx2);
J(ip,ipyp)=J(ip,ipyp)+2/dy2/(dy1+dy2);
J(ip,ipym)=J(ip,ipym)+2/dy1/(dy1+dy2);
    end
end

%corners
ix=1; iy=1; ip=iy+(ix-1)*ny; g(ip)=u(ip); J(ip,ip)=1;
ix=1; iy=ny; ip=iy+(ix-1)*ny; g(ip)=u(ip); J(ip,ip)=1;
ix=nx; iy=1; ip=iy+(ix-1)*ny; g(ip)=u(ip); J(ip,ip)=1;
ix=nx; iy=ny; ip=iy+(ix-1)*ny; g(ip)=u(ip); J(ip,ip)=1;

%bords
    for iy=1+1:length(y)-1
        ix=1; ip=iy+(ix-1)*ny; g(ip)=u(ip); J(ip,ip)=1;
        ix=nx; ip=iy+(ix-1)*ny; g(ip)=u(ip); J(ip,ip)=1;
    end
    for ix=1+1:length(x)-1
        iy=1; ip=iy+(ix-1)*ny; g(ip)=u(ip); J(ip,ip)=1;
        iy=ny; ip=iy+(ix-1)*ny; g(ip)=u(ip); J(ip,ip)=1;
    end

%
u=u-J\g;
%
uPhys=reshape(u,[ny,nx]);
mesh(x,y,uPhys);
%
ix=(length(x)+1)/2; iy=(length(y)+1)/2; ip=iy+(ix-1)*ny; 
uarr=[uarr;u(ip)]
% anal=-1/pi/pi*2
%%
d=diff(uarr)
d(1)/d(2)
