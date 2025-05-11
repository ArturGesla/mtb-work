clc; clear;
uarr=[];
%
nx=5; 
narr=[];
%%
% nx=(nx-1)*2+1;
ny=nx;
x=linspace(-1,1,nx)'; y=x;
dx=mean(diff(x));

%d2udx2-1=0
g=zeros(nx*ny*3,1); %laplacian d2x d2y
u=g-1; %u d2udx2 d2udy2
J=sparse(nx*ny*3,nx*ny*3,0);
%
% f=x*0-1*x.^2+1;
% f=x*0-1*x+1;
f=u*0-1;
% f=cos(pi*x);
%%

g=zeros(nx*ny*3,1);
% u=g-1;
J=sparse(nx*ny*3,nx*ny*3,0);

ii=zeros(4*nx*ny,1); jj=ii; vv=ii; iip=1;

for ix=1+1:length(x)-1
    for iy=1+1:length(y)-1
        
%laplacian
ip=(iy+(ix-1)*ny-1)*3+1;
g(ip)=g(ip)+u(ip+1)+u(ip+2)-f(ip)*cos(pi/2*x(ix))*cos(pi/2*y(iy));
ii(iip)=ip; jj(iip)=ip+1; vv(iip)=1; iip=iip+1;
ii(iip)=ip; jj(iip)=ip+2; vv(iip)=1; iip=iip+1;

%implicit eq for second x der
ip=(iy+(ix-1)*ny-1)*3+2;
ipp=(iy+(ix+1-1)*ny-1)*3+2;
ipm=(iy+(ix-1-1)*ny-1)*3+2;
ipf=(iy+(ix-1)*ny-1)*3+1;
ipfp=(iy+(ix+1-1)*ny-1)*3+1;
ipfm=(iy+(ix-1-1)*ny-1)*3+1;
g(ip)=g(ip)+1/10*u(ipp)+u(ip)+1/10*u(ipm)-6/5/dx/dx*(u(ipfp)-2*u(ipf)+u(ipfm));
ii(iip)=ip; jj(iip)=ipp; vv(iip)=1/10; iip=iip+1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;
ii(iip)=ip; jj(iip)=ipm; vv(iip)=1/10; iip=iip+1;
ii(iip)=ip; jj(iip)=ipfp; vv(iip)=-6/5/dx/dx; iip=iip+1;
ii(iip)=ip; jj(iip)=ipf; vv(iip)=-6/5/dx/dx*-2; iip=iip+1;
ii(iip)=ip; jj(iip)=ipfm; vv(iip)=-6/5/dx/dx; iip=iip+1;

%implicit eq for second y der
ip=(iy+(ix-1)*ny-1)*3+3;
ipp=(iy+1+(ix-1)*ny-1)*3+3;
ipm=(iy-1+(ix-1)*ny-1)*3+3;
ipf=(iy+(ix-1)*ny-1)*3+1;
ipfp=(iy+1+(ix-1)*ny-1)*3+1;
ipfm=(iy-1+(ix-1)*ny-1)*3+1;
g(ip)=g(ip)+1/10*u(ipp)+u(ip)+1/10*u(ipm)-6/5/dx/dx*(u(ipfp)-2*u(ipf)+u(ipfm));
ii(iip)=ip; jj(iip)=ipp; vv(iip)=1/10; iip=iip+1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;
ii(iip)=ip; jj(iip)=ipm; vv(iip)=1/10; iip=iip+1;
ii(iip)=ip; jj(iip)=ipfp; vv(iip)=-6/5/dx/dx; iip=iip+1;
ii(iip)=ip; jj(iip)=ipf; vv(iip)=-6/5/dx/dx*-2; iip=iip+1;
ii(iip)=ip; jj(iip)=ipfm; vv(iip)=-6/5/dx/dx; iip=iip+1;






% 
% 
% 
% 
% ipxp=iy+(ix-1+1)*ny;
% ipxm=iy+(ix-1-1)*ny;
% ipyp=iy+(ix-1)*ny+1;
% ipym=iy+(ix-1)*ny-1;
% 
% g(ip)=g(ip)+-2/dx1/dx2*u(ip)+2/dx2/(dx1+dx2)*u(ipxp)+2/dx1/(dx1+dx2)*u(ipxm);
% g(ip)=g(ip)+-2/dy1/dy2*u(ip)+2/dy2/(dy1+dy2)*u(ipyp)+2/dy1/(dy1+dy2)*u(ipym);
% g(ip)=g(ip)+-f(ip)*cos(pi/2*x(ix))*cos(pi/2*y(iy));
% % g(ip)=g(ip)+-f(ip);
% % g(ip)=g(ip)+-u(ip)*pi*pi/4;
% 
% % J(ip,ip)=J(ip,ip)+-2/dx1/dx2-2/dy1/dy2;
% % % J(ip,ip)=J(ip,ip)+-1*pi*pi/4;
% % J(ip,ipxp)=J(ip,ipxp)+2/dx2/(dx1+dx2);
% % J(ip,ipxm)=J(ip,ipxm)+2/dx1/(dx1+dx2);
% % J(ip,ipyp)=J(ip,ipyp)+2/dy2/(dy1+dy2);
% % J(ip,ipym)=J(ip,ipym)+2/dy1/(dy1+dy2);
% 
% 
% ii(iip)=ip; jj(iip)=ip; vv(iip)=-2/dx1/dx2-2/dy1/dy2; iip=iip+1;
% ii(iip)=ip; jj(iip)=ipxp; vv(iip)=+2/dx2/(dx1+dx2); iip=iip+1;
% ii(iip)=ip; jj(iip)=ipxm; vv(iip)=+2/dx1/(dx1+dx2); iip=iip+1;
% ii(iip)=ip; jj(iip)=ipyp; vv(iip)=+2/dy2/(dy1+dy2); iip=iip+1;
% ii(iip)=ip; jj(iip)=ipym; vv(iip)=+2/dy1/(dy1+dy2); iip=iip+1;


    end
end

%%

%corners
ix=1; iy=1; 
ip=(iy+(ix-1)*ny-1)*3+1; %NOT ENOUGH EQS!!
% ip=iy+(ix-1)*ny; 
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
        ix=1; ip=iy+(ix-1)*ny; g(ip)=u(ip); 
%         J(ip,ip)=1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;

        ix=nx; ip=iy+(ix-1)*ny; g(ip)=u(ip); 
%         J(ip,ip)=1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;

    end
    for ix=1+1:length(x)-1
        iy=1; ip=iy+(ix-1)*ny; g(ip)=u(ip); 
%         J(ip,ip)=1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;

        iy=ny; ip=iy+(ix-1)*ny; g(ip)=u(ip);
%         J(ip,ip)=1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;

    end

    J=J+sparse([ii(1:iip-1);nx*ny],[jj(1:iip-1);nx*ny],[vv(1:iip-1);0]);


%
tic;
u=u-J\g;
toc;
norm(g)
%
%
uPhys=reshape(u,[ny,nx]);
mesh(x,y,uPhys);
%
ix=(length(x)+1)/2; iy=(length(y)+1)/2; ip=iy+(ix-1)*ny; 
% uarr=[uarr;u(ip)]
% narr=[narr;nx];
% nx=(nx-1)*2+1;
% anal=-1/pi/pi*2
%
% clf; 
B=zeros(nx,ny); B(2:end-1,2:end-1)=1; B=reshape(B,[nx*ny,1]); B=sparse(1:length(B),1:length(B),B);
% [ev,evs]=eig(full(J)); evs=diag(evs); 
[ev,evs]=eigs((J),B,20,2); evs=diag(evs);
save('evsFD','evs')
% plot(evs,'x');

uarr=[uarr;real(evs(1:5)')]
narr=[narr;nx];
%%
save('evsFD',"uarr","narr");

%%
% i=1; 
u=real(ev(:,i));
uPhys=reshape(u,[ny,nx]);
mesh(x,y,uPhys);
i=i+1;
%%
d=diff(uarr)
d(1)/d(2)
%%
clf;
loglog(narr,abs(uarr-2/pi/pi),'x-'); hold on;
loglog(narr,narr.^(-2),'-')
