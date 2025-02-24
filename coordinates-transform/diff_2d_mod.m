clc; clear;
uarr=[];
%
% nx=71;
nx=31;
% nx=71*2+1;
% ny=(nx+1)/2;
ny=nx-10;
narr=[];
L=2;
%
% nx=(nx-1)*2+1;
% ny=nx;
% x=linspace(0,1,nx)';
% y=linspace(0,1,ny)';

y=linspace(0,1,ny); x=linspace(0,atan(0.5),nx); 
%
dx=mean(diff(x));
dy=mean(diff(y));

%d2udx2-1=0
g=zeros(nx*ny,1); %laplacian d2x+d2y
u=g; %u d2udx2+d2udy2
J=sparse(nx*ny,nx*ny,0);
%
% f=x*0-1*x.^2+1;
% f=x*0-1*x+1;
f=u*0;
% f=cos(pi*x);
%
%
g=zeros(nx*ny,1);
% u=g-1;
J=sparse(nx*ny,nx*ny,0);
B=J;

ii=zeros(4*nx*ny,1); jj=ii; vv=ii; iip=1;
iib=zeros(4*nx*ny,1); jjb=iib; vvb=iib; iipb=1;
%
for ix=1+1:length(x)-1
    for iy=1+1:length(y)-1
            % %second order
            ipf=(iy+(ix-1)*ny-1)+1;
            ipfxp=(iy+(ix+1-1)*ny-1)+1;
            ipfxm=(iy+(ix-1-1)*ny-1)+1;
            ipfyp=(iy+1+(ix-1)*ny-1)+1;
            ipfym=(iy-1+(ix-1)*ny-1)+1;
            ip=ipf;

            th=x(ix);
            z=y(iy);
            val1=cos(th)^2/(L-z)*cos(th)^2/(L-z);
            val2=cos(th)^2/(L-z)*(-2*cos(th)*sin(th)/(L-z));
%             val2=cos(th)^2/(L-z)*(-10*cos(th)*sin(th)/(L-z))*0;

            g(ip)=g(ip)+(u(ipfxp)+u(ipfxm)-2*u(ipf))/dx/dx*val1;
            g(ip)=g(ip)+(u(ipfxp)-u(ipfxm))/dx/2*val2;
            g(ip)=g(ip)+(u(ipfyp)+u(ipfym)-2*u(ipf))/dy/dy;

            ii(iip)=ip; jj(iip)=ipfxp; vv(iip)=1/dx/dx*val1; iip=iip+1;
            ii(iip)=ip; jj(iip)=ipfxm; vv(iip)=1/dx/dx*val1; iip=iip+1;
                        ii(iip)=ip; jj(iip)=ipf; vv(iip)=-2/dx/dx*val1; iip=iip+1;


            ii(iip)=ip; jj(iip)=ipfxp; vv(iip)=1/dx/2*val2; iip=iip+1;
            ii(iip)=ip; jj(iip)=ipfxm; vv(iip)=-1/dx/2*val2; iip=iip+1;

            ii(iip)=ip; jj(iip)=ipfyp; vv(iip)=1/dy/dy; iip=iip+1;
            ii(iip)=ip; jj(iip)=ipfym; vv(iip)=1/dy/dy; iip=iip+1;
            ii(iip)=ip; jj(iip)=ipf; vv(iip)=-2/dy/dy; iip=iip+1;

%             u(ip)=2;
        % iib(iipb)=ipf; jjb(iipb)=ipf; vvb(iipb)=1; iipb=iipb+1;

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

% %forcing
% for ix=1:length(x)
%     for iy=1:length(y)
%         ipf=(iy+(ix-1)*ny-1)+1;
% f(ipf)=-cos(pi/2*x(ix))*cos(pi/2*y(iy));
%     end
% end
%
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
    %left
    ix=1; ip=iy+(ix-1)*ny; g(ip)=u(ip);
    %         J(ip,ip)=1;
    ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;

    ix=nx; ip=iy+(ix-1)*ny; g(ip)=u(ip)-1;
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
%
J=J+sparse([ii(1:iip-1);nx*ny],[jj(1:iip-1);nx*ny],[vv(1:iip-1);0]);
% B=B+sparse([iib(1:iipb-1);nx*ny],[jjb(1:iipb-1);nx*ny],[vvb(1:iipb-1);0]);
%

%
tic;
u=u-J\g;
% u=J\(B*f);
toc;
% norm(g)
%
% plt
uPhys=reshape(u,[ny,nx]);
% mesh(x,y,uPhys);
mesh(1:nx,1:ny,uPhys);
%
% pcolor(1:nx,1:ny,uPhys); shading interp; colormap(parula(8))
% contour(1:nx,1:ny,uPhys,30); 
%
%
[TH,Z]=meshgrid(x,y);
X=(2-Z).*tan(TH); Y=Z;

%
clf;
% mesh(X,Y,uPhys)
contour(X,Y,uPhys,[0:0.1:1],'r'); 

%
clf; 
a=load("dataCart.mat");
a.nx
[X2,Y2]=meshgrid(a.x,a.y);
mesh(X2,Y2,a.uPhys,'EdgeColor','r');
hold on;
mesh(X,Y,uPhys,'EdgeColor','b')

%
pp=interp2(X2,Y2,a.uPhys,X,Y);
max(max(uPhys-pp))
%%

clf;
contour(X,Y,uPhys,[0:0.1:1],'r');
hold on;
contour(X2,Y2,a.uPhys,[0:0.1:1],'k');
legend("modified","cartesian"); grid on;
title("contour 0:0.1:1 | nth or nx: "+num2str(nx));
% exportgraphics(gcf,"comparison-"+num2str(nx)+".png")
%nie wiem czemu sie nie zgadza 24/02/2025
%%
clf;
mesh(X,Y,uPhys-pp,'EdgeColor','b'); hold on;
% mesh(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),uPhys(2:end-1,2:end-1)-pp(2:end-1,2:end-1),'EdgeColor','b'); hold on;
% mesh(X,Y,pp,'EdgeColor','r')


