clc; clear;
uarr=[];
%%
nt=21;
nx=nt; ny=nx;
% x=linspace(-1,1,nx)';
% x=x.^3;
x=-cos(linspace(0,pi,nx))'; y=x;
%d2udx2-1=0
%
g=zeros(nx*ny,1);
u=g;
J=sparse(nx*ny,nx*ny,0);
%
% f=x*0-1*x.^2+1;
% f=x*0-1*x+1;
f=u*0-1;

%
for ix=1+1:length(x)-1
    for iy=1+1:length(y)-1
        ip=iy+(ix-1)*ny; xc=x(ix); yc=y(iy);

        for ikx=0:nt-1
            for iky=0:nt-1
                ipk=iky+1+(ikx-1+1)*ny;
                g(ip)=g(ip)+u(ipk)*(dTn2(xc,ikx)*Tn(yc,iky)+dTn2(yc,iky)*Tn(xc,ikx));
                J(ip,ipk)=J(ip,ipk)+(dTn2(xc,ikx)*Tn(yc,iky)+dTn2(yc,iky)*Tn(xc,ikx));
            end
        end
        g(ip)=g(ip)-f(ip);


    end
end
%
%corners
ix=1; iy=1; ip=iy+(ix-1)*ny;
for ikx=0:nt-1
    for iky=0:nt-1
        ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
        g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
        J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
    end
end
ix=1; iy=ny; ip=iy+(ix-1)*ny;
for ikx=0:nt-1
    for iky=0:nt-1
        ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
        g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
        J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
    end
end
ix=nx; iy=1; ip=iy+(ix-1)*ny;
for ikx=0:nt-1
    for iky=0:nt-1
        ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
        g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
        J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
    end
end
ix=nx; iy=ny; ip=iy+(ix-1)*ny;
for ikx=0:nt-1
    for iky=0:nt-1
        ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
        g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
        J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
    end
end
%
%bords
for iy=1+1:length(y)-1
    ix=1; ip=iy+(ix-1)*ny;
    for ikx=0:nt-1
        for iky=0:nt-1
            ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
            g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
            J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
        end
    end
    ix=nx; ip=iy+(ix-1)*ny;
    for ikx=0:nt-1
        for iky=0:nt-1
            ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
            g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
            J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
        end
    end
end
for ix=1+1:length(x)-1
    iy=1; ip=iy+(ix-1)*ny;
    for ikx=0:nt-1
        for iky=0:nt-1
            ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
            g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
            J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
        end
    end
    iy=ny; ip=iy+(ix-1)*ny;
    for ikx=0:nt-1
        for iky=0:nt-1
            ipk=iky+1+(ikx-1+1)*ny;xc=x(ix); yc=y(iy);
            g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
            J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
        end
    end
end
%
u=u-J\g;
%
uPhys=zeros(length(x)*length(y),1);
for ix=1:length(x)
    for iy=1:length(y)
        ip=iy+(ix-1)*length(y); xc=x(ix); yc=y(iy);

        for ikx=0:nt-1
            for iky=0:nt-1
                ipk=iky+1+(ikx-1+1)*ny;
                uPhys(ip)=uPhys(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
            end
        end
    end
end
%
ix=(length(x)+1)/2; iy=(length(y)+1)/2; ip=iy+(ix-1)*ny; 
uarr=[uarr;uPhys(ip)]
%%
mesh(reshape(log10(abs(u)),[ny,nx]))
%%
d=diff(uarr)
d(1)/d(2)
%%
% x=-cos(linspace(0,pi,nx*20))'; y=x;

% uarr=[uarr;u((length(x)+1)/2)]
uPhys=zeros(length(x)*length(y),1);
for ix=1:length(x)
    for iy=1:length(y)
        ip=iy+(ix-1)*length(y); xc=x(ix); yc=y(iy);

        for ikx=0:nt-1
            for iky=0:nt-1
                ipk=iky+1+(ikx-1+1)*ny;
                uPhys(ip)=uPhys(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
            end
        end
    end
end
%%
uPhys2=reshape(uPhys,[length(y),length(x)]);
mesh(x,y,uPhys2);
%%
uarr=[uarr;uPhys((length(x)+1)/2)]
%%
d=diff(uarr)
d(1)/d(2)
%%
close all;
plot(x,uPhys,'-x')