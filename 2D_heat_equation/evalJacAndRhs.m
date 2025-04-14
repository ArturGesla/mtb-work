function [g,jac]=evalJacAndRhs(x,u,f,k)


y=x;
nx=length(x);
ny=length(y);

g=zeros(nx*ny,1); %laplacian d2x+d2y
% u=g; %u d2udx2+d2udy2
J=sparse(nx*ny,nx*ny,0);

dx=mean(diff(x));

B=J;



ii=zeros(4*nx*ny,1); jj=ii; vv=ii; iip=1;
iib=zeros(4*nx*ny,1); jjb=iib; vvb=iib; iipb=1;

for ix=1+1:length(x)-1
    for iy=1+1:length(y)-1
        
% %second order
ipf=(iy+(ix-1)*ny-1)+1;
ipfxp=(iy+(ix+1-1)*ny-1)+1;
ipfxm=(iy+(ix-1-1)*ny-1)+1;
ipfyp=(iy+1+(ix-1)*ny-1)+1;
ipfym=(iy-1+(ix-1)*ny-1)+1;
ip=ipf;

g(ip)=g(ip)+(u(ipfxp)+u(ipfxm)+u(ipfyp)+u(ipfym)-4*u(ipf))/dx/dx;
ii(iip)=ip; jj(iip)=ipfxp; vv(iip)=1/dx/dx; iip=iip+1;
ii(iip)=ip; jj(iip)=ipfxm; vv(iip)=1/dx/dx; iip=iip+1;
ii(iip)=ip; jj(iip)=ipfyp; vv(iip)=1/dx/dx; iip=iip+1;
ii(iip)=ip; jj(iip)=ipfym; vv(iip)=1/dx/dx; iip=iip+1;
ii(iip)=ip; jj(iip)=ipf; vv(iip)=-4/dx/dx; iip=iip+1;

iib(iipb)=ipf; jjb(iipb)=ipf; vvb(iipb)=1; iipb=iipb+1;

%forcing
% g(ip)=g(ip)+rand();
g(ip)=g(ip)+f(ip);

%polar
g(ip)=g(ip)+1/x(ix)*(u(ipfxp)-u(ipfxm))/2/dx;
ii(iip)=ip; jj(iip)=ipfxp; vv(iip)=1/dx/2/x(ix); iip=iip+1;
ii(iip)=ip; jj(iip)=ipfxm; vv(iip)=-1/dx/2/x(ix); iip=iip+1;

%polar k
g(ip)=g(ip)-1/x(ix)/x(ix)*k*k*(u(ipf));
ii(iip)=ip; jj(iip)=ipf; vv(iip)=-1/x(ix)/x(ix)*k*k; iip=iip+1;

%fourth order what a magic
% ipf=(iy+(ix-1)*ny-1)+1;
% ipfxp=(iy+(ix+1-1)*ny-1)+1;
% ipfxpyp=(iy+1+(ix+1-1)*ny-1)+1;
% ipfxpym=(iy-1+(ix+1-1)*ny-1)+1;
% ipfxmyp=(iy+1+(ix-1-1)*ny-1)+1;
% ipfxmym=(iy-1+(ix-1-1)*ny-1)+1;
% ipfxm=(iy+(ix-1-1)*ny-1)+1;
% ipfyp=(iy+1+(ix-1)*ny-1)+1;
% ipfym=(iy-1+(ix-1)*ny-1)+1;
% ip=ipf;
% 
% % g(ip)=g(ip)+(u(ipfxp)+u(ipfxm)+u(ipfyp)+u(ipfym)-4*u(ipf))/dx/dx;
% ii(iip)=ip; jj(iip)=ipfxp; vv(iip)=4/6/dx/dx; iip=iip+1;
% ii(iip)=ip; jj(iip)=ipfxpyp; vv(iip)=1/6/dx/dx; iip=iip+1;
% ii(iip)=ip; jj(iip)=ipfxpym; vv(iip)=1/6/dx/dx; iip=iip+1;
% ii(iip)=ip; jj(iip)=ipfxm; vv(iip)=4/6/dx/dx; iip=iip+1;
% ii(iip)=ip; jj(iip)=ipfxmyp; vv(iip)=1/6/dx/dx; iip=iip+1;
% ii(iip)=ip; jj(iip)=ipfxmym; vv(iip)=1/6/dx/dx; iip=iip+1;
% ii(iip)=ip; jj(iip)=ipfyp; vv(iip)=4/6/dx/dx; iip=iip+1;
% ii(iip)=ip; jj(iip)=ipfym; vv(iip)=4/6/dx/dx; iip=iip+1;
% ii(iip)=ip; jj(iip)=ipf; vv(iip)=-20/6/dx/dx; iip=iip+1;
% 
% 
% 
% iib(iipb)=ipf; jjb(iipb)=ipf; vvb(iipb)=8/12; iipb=iipb+1;
% iib(iipb)=ipf; jjb(iipb)=ipfxp; vvb(iipb)=1/12; iipb=iipb+1;
% iib(iipb)=ipf; jjb(iipb)=ipfxm; vvb(iipb)=1/12; iipb=iipb+1;
% iib(iipb)=ipf; jjb(iipb)=ipfyp; vvb(iipb)=1/12; iipb=iipb+1;
% iib(iipb)=ipf; jjb(iipb)=ipfym; vvb(iipb)=1/12; iipb=iipb+1;


    end
end

%forcing
for ix=1:length(x)
    for iy=1:length(y)
        ipf=(iy+(ix-1)*ny-1)+1;
f(ipf)=-cos(pi/2*x(ix))*cos(pi/2*y(iy));
    end
end
%

%corners
ix=1; iy=1; ip=iy+(ix-1)*ny; 
g(ip)=u(ip); 
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;

ix=1; iy=ny; ip=iy+(ix-1)*ny; g(ip)=u(ip);
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;

ix=nx; iy=1; ip=iy+(ix-1)*ny; g(ip)=u(ip); 
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;

ix=nx; iy=ny; ip=iy+(ix-1)*ny; g(ip)=u(ip);
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;


%bords
    for iy=1+1:length(y)-1
        ix=1; ip=iy+(ix-1)*ny; g(ip)=u(ip)-u(ip+ny); 
%         ix=1; ip=iy+(ix-1)*ny; g(ip)=u(ip)+u(ip+ny); 
%         J(ip,ip)=1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;
% ii(iip)=ip; jj(iip)=ip+ny; vv(iip)=1; iip=iip+1;
ii(iip)=ip; jj(iip)=ip+ny; vv(iip)=-1; iip=iip+1;

%         ix=nx; ip=iy+(ix-1)*ny; g(ip)=u(ip)+u(ip-ny)-1; 
%         ix=nx; ip=iy+(ix-1)*ny; g(ip)=u(ip)+u(ip-ny); 
        ix=nx; ip=iy+(ix-1)*ny; g(ip)=u(ip)-u(ip-ny); 
%         J(ip,ip)=1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;
% ii(iip)=ip; jj(iip)=ip-ny; vv(iip)=1; iip=iip+1;
ii(iip)=ip; jj(iip)=ip-ny; vv(iip)=-1; iip=iip+1;

    end
    for ix=1+1:length(x)-1
%         iy=1; ip=iy+(ix-1)*ny; g(ip)=u(ip)+u(ip+1); 
        iy=1; ip=iy+(ix-1)*ny; g(ip)=u(ip)-u(ip+1); 
%         J(ip,ip)=1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;
% ii(iip)=ip; jj(iip)=ip+1; vv(iip)=1; iip=iip+1;
ii(iip)=ip; jj(iip)=ip+1; vv(iip)=-1; iip=iip+1;

%         iy=ny; ip=iy+(ix-1)*ny; g(ip)=u(ip)+u(ip-1);
        iy=ny; ip=iy+(ix-1)*ny; g(ip)=u(ip)-u(ip-1);
%         J(ip,ip)=1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;
% ii(iip)=ip; jj(iip)=ip-1; vv(iip)=1; iip=iip+1;
ii(iip)=ip; jj(iip)=ip-1; vv(iip)=-1; iip=iip+1;

    end

    J=J+sparse([ii(1:iip-1);nx*ny],[jj(1:iip-1);nx*ny],[vv(1:iip-1);0]);
    B=B+sparse([iib(1:iipb-1);nx*ny],[jjb(1:iipb-1);nx*ny],[vvb(1:iipb-1);0]);

jac=J;
end