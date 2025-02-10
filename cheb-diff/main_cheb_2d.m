clc; clear;
uarr=[];
narr=[];
%
nt=3;
%%
nt=nt+2;
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
%         g(ip)=g(ip)-f(ip);
%         g(ip)=g(ip)-f(ip)*cos(pi/2*xc)*cos(pi/2*yc);
%         g(ip)=g(ip)-f(ip)*cos(pi*xc)*cos(pi/2*yc);


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
%     g(ip)=g(ip)-1; % bc gives 2nd order maybe
%     g(ip)=g(ip)-(1-y(iy).^2); % bc maybe some 5th
%     g(ip)=g(ip)-(1-abs(y(iy))); % bc maybe 3rd
    g(ip)=g(ip)-exp(-1./(1-y(iy).^2)); % bc

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

%     g(ip)=g(ip)-1; % bc

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

uecc=0;
        for ikx=0:nt-1
            for iky=0:nt-1
                ipk=iky+1+(ikx-1+1)*ny;
                uecc=uecc+u(ipk)*(Tn(0.6,iky)*Tn(0.6,ikx));
            end
        end


%
ix=(length(x)+1)/2; iy=(length(y)+1)/2; ip=iy+(ix-1)*ny; 
% uarr=[uarr;uPhys(ip)]
uarr=[uarr;uecc]
narr=[narr; nt];

%%
mesh(reshape(log10(abs(u)),[ny,nx]))
% mesh(reshape((abs(u)),[ny,nx]))
%%

clf;
% semilogy(narr,abs(uarr-2/pi/pi),'x-'); hold on;
% semilogy(narr,abs(uarr-uarr(end)),'x-'); hold on;
loglog(narr,abs(uarr-uarr(end)),'x-'); hold on;
% semilogy(narr,narr.^(-2),'o-')
loglog(narr,narr.^(-2),'o-')
loglog(narr,narr.^(-3),'o-')
loglog(narr,narr.^(-4),'o-')

%%
d=diff(uarr)
d(1)/d(2)
%%
x=-cos(linspace(0,pi,nx*5))'; y=x;
% x=-cos(linspace(0,pi,nx))'; y=x;

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
%
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