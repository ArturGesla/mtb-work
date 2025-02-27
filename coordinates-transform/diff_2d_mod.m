clc; clear;
uarr=[];
%
% nx=71;
% nx=301;
% nx=71*2+1;
% ny=(nx+1)/2;
nx=160+1;
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


        ipfxpyp=(iy+1+(ix+1-1)*ny-1)+1;
        ipfxmym=(iy-1+(ix-1-1)*ny-1)+1;
        ipfxmyp=(iy+1+(ix-1-1)*ny-1)+1;
        ipfxpym=(iy-1+(ix-1+1)*ny-1)+1;

        ip=ipf;

        th=x(ix);
        z=y(iy);

        % d2fdx2
        val1=cos(th)^2/(L-z)*cos(th)^2/(L-z);
        val2=cos(th)^2/(L-z)*(-2*cos(th)*sin(th)/(L-z));
        %
        g(ip)=g(ip)+(u(ipfxp)+u(ipfxm)-2*u(ipf))/dx/dx*val1;
        g(ip)=g(ip)+(u(ipfxp)-u(ipfxm))/dx/2*val2;
%         g(ip)=g(ip)+(u(ipfyp)+u(ipfym)-2*u(ipf))/dy/dy;

        ii(iip)=ip; jj(iip)=ipfxp; vv(iip)=1/dx/dx*val1; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipfxm; vv(iip)=1/dx/dx*val1; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipf; vv(iip)=-2/dx/dx*val1; iip=iip+1;


        ii(iip)=ip; jj(iip)=ipfxp; vv(iip)=1/dx/2*val2; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipfxm; vv(iip)=-1/dx/2*val2; iip=iip+1;

%         ii(iip)=ip; jj(iip)=ipfyp; vv(iip)=1/dy/dy; iip=iip+1;
%         ii(iip)=ip; jj(iip)=ipfym; vv(iip)=1/dy/dy; iip=iip+1;
%         ii(iip)=ip; jj(iip)=ipf; vv(iip)=-2/dy/dy; iip=iip+1;

%d2fdy2
        val3=sin(th)*cos(th)/(L-z)^2;
        val4=sin(th)*cos(th)/(L-z)*2;
        val5=val4*cos(2*th)/(L-z)/2;
        val6=val4^2/2/2;
        
        %
        g(ip)=g(ip)+(u(ipfyp)+u(ipfym)-2*u(ipf))/dy/dy;
        g(ip)=g(ip)+(u(ipfxp)-u(ipfxm))/dx/2*val3;
        g(ip)=g(ip)+((u(ipfxpyp)-u(ipfxpym))/dy/2-(u(ipfxmyp)-u(ipfxmym))/dy/2)/2/dx*val4;
        g(ip)=g(ip)+(u(ipfxp)-u(ipfxm))/dx/2*val5;
        g(ip)=g(ip)+(u(ipfxp)+u(ipfxm)-2*u(ipf))/dx/dx*val6;
        

        ii(iip)=ip; jj(iip)=ipfyp; vv(iip)=1/dy/dy; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipfym; vv(iip)=1/dy/dy; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipf; vv(iip)=-2/dy/dy; iip=iip+1;

        ii(iip)=ip; jj(iip)=ipfxp; vv(iip)=1/dx/2*val3; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipfxm; vv(iip)=-1/dx/2*val3; iip=iip+1;

        ii(iip)=ip; jj(iip)=ipfxpyp; vv(iip)=1/dx/2/2/dy*val4; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipfxpym; vv(iip)=-1/dx/2/2/dy*val4; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipfxmyp; vv(iip)=-1/dx/2/2/dy*val4; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipfxmym; vv(iip)=1/dx/2/2/dy*val4; iip=iip+1;

        ii(iip)=ip; jj(iip)=ipfxp; vv(iip)=1/dx/2*val5; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipfxm; vv(iip)=-1/dx/2*val5; iip=iip+1;

        ii(iip)=ip; jj(iip)=ipfxp; vv(iip)=1/dx/dx*val6; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipfxm; vv(iip)=1/dx/dx*val6; iip=iip+1;
        ii(iip)=ip; jj(iip)=ipf; vv(iip)=-2/dx/dx*val6; iip=iip+1;


       
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

ix=nx; iy=1; ip=iy+(ix-1)*ny; g(ip)=u(ip)-1;
ii(iip)=ip; jj(iip)=ip; vv(iip)=1; iip=iip+1;

ix=nx; iy=ny; ip=iy+(ix-1)*ny; g(ip)=u(ip)-1;
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
uPhys=reshape(u,[ny,nx]);

[TH,Z]=meshgrid(x,y);
X=(2-Z).*tan(TH); Y=Z;
save("dataMod-"+num2str(nx),'X','Y','uPhys','nx');

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

%
clf;
% mesh(X,Y,uPhys)
contour(X,Y,uPhys,[0:0.1:1],'r');

%
clf;
a=load("dataCart-81.mat");
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
%%
clf;
% mesh(X,Y,uPhys-pp,'EdgeColor','b'); hold on;
% mesh(X(2:end-1,2:end-1),Y(2:end-1,2:end-1),uPhys(2:end-1,2:end-1)-pp(2:end-1,2:end-1),'EdgeColor','b'); hold on;
% mesh(X(3:end-2,3:end-2),Y(3:end-2,3:end-2),uPhys(3:end-2,3:end-2)-pp(3:end-2,3:end-2),'EdgeColor','b'); hold on;
mesh(X(3:end-2,3:end-2),Y(3:end-2,3:end-2),uPhys(3:end-2,3:end-2)-pp(3:end-2,3:end-2)); hold on;
% mesh(X,Y,pp,'EdgeColor','r')
% mesh(X,Y,uPhys,'EdgeColor','r')


%%  conv  cart


clf; hold on;
% a=load("dataCart-21.mat"); contour(a.x,a.y,a.uPhys,[0:0.1:1],'k');
% a=load("dataCart-41.mat"); contour(a.x,a.y,a.uPhys,[0:0.1:1],'r');
% a=load("dataCart-81.mat"); contour(a.x,a.y,a.uPhys,[0:0.1:1],'b');
% a=load("dataCart-161.mat"); contour(a.x,a.y,a.uPhys,[0:0.1:1],'r');
a=load("dataCart-321.mat"); contour(a.x,a.y,a.uPhys,[0:0.1:1],'k');

% a=load("dataMod-21.mat"); contour(a.X,a.Y,a.uPhys,[0:0.1:1],'k');
% a=load("dataMod-41.mat"); contour(a.X,a.Y,a.uPhys,[0:0.1:1],'r');
% a=load("dataMod-81.mat"); contour(a.X,a.Y,a.uPhys,[0:0.1:1],'b');
a=load("dataMod-161.mat"); contour(a.X,a.Y,a.uPhys,[0:0.1:1],'r');
%works nice!

