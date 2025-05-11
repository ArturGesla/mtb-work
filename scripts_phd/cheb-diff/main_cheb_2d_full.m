clc; clear;
uarr=[];
narr=[];
%
nt=3;
%
% for int=1:20
% for nt=[3:2:41,51]
for nt=71%[3:2:41,51]
% nt=nt+2;
nx=nt; ny=nx;
% x=linspace(-1,1,nx)';
% x=x.^3;
x=-cos(linspace(0,pi,nx))'; y=x;
[X,Y]=meshgrid(x,y);
%d2udx2-1=0
%
g=zeros(nx*ny,1);
u=g;
J=zeros(nx*ny,nx*ny);
B=zeros(nx*ny,nx*ny);
%
% f=x*0-1*x.^2+1;
% f=x*0-1*x+1;
f=u*0-1;
%
tic; 
% ii=zeros(nt*nt*nx*ny,1); jj=ii; vv=ii; iip=1;
% iib=zeros(nt*nt*nx*ny,1); jjb=iib; vvb=iib; iipb=1;
for ix=1+1:length(x)-1
    for iy=1+1:length(y)-1
        ip=iy+(ix-1)*ny; xc=x(ix); yc=y(iy);

        for ikx=0:nt-1
            for iky=0:nt-1
                ipk=iky+1+(ikx-1+1)*ny;
                aa=(dTn2(xc,ikx)*Tn(yc,iky)+dTn2(yc,iky)*Tn(xc,ikx));
%                 g(ip)=g(ip)+u(ipk)*(dTn2(xc,ikx)*Tn(yc,iky)+dTn2(yc,iky)*Tn(xc,ikx));
                g(ip)=g(ip)+u(ipk)*aa;
                J(ip,ipk)=J(ip,ipk)+(dTn2(xc,ikx)*Tn(yc,iky)+dTn2(yc,iky)*Tn(xc,ikx));
% ii(iip)=ip; jj(iip)=ipk; vv(iip)=(dTn2(xc,ikx)*Tn(yc,iky)+dTn2(yc,iky)*Tn(xc,ikx)); iip=iip+1;
% ii(iip)=ip; jj(iip)=ipk; vv(iip)=aa; iip=iip+1;
% iib(iipb)=ip; jjb(iipb)=ipk; vvb(iipb)=Tn(xc,ikx)*Tn(yc,iky); iipb=iipb+1;
            end
        end
        %focring
%         g(ip)=g(ip)-f(ip);
%         g(ip)=g(ip)-f(ip)*cos(pi/2*xc)*cos(pi/2*yc);
%         g(ip)=g(ip)-f(ip)*cos(pi*xc)*cos(pi/2*yc);


    end
end

% J=J+sparse([ii(1:iip-1);nx*ny],[jj(1:iip-1);nx*ny],[vv(1:iip-1);0]);
% B=B+sparse([iib(1:iipb-1);nx*ny],[jjb(1:iipb-1);nx*ny],[vvb(1:iipb-1);0]);
toc;
% gold=g;

%
% tic;
%
%corners
ii=zeros(nt*nt*nx*ny,1); jj=ii; vv=ii; iip=1;

ix=1; iy=1; ip=iy+(ix-1)*ny;
for ikx=0:nt-1
    for iky=0:nt-1
        ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
        g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
%         J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
        ii(iip)=ip; jj(iip)=ipk; vv(iip)=(Tn(yc,iky)*Tn(xc,ikx)); iip=iip+1;

    end
end
% g(ip)=g(ip)-1;

ix=1; iy=ny; ip=iy+(ix-1)*ny;
for ikx=0:nt-1
    for iky=0:nt-1
        ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
        g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
%         J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
        ii(iip)=ip; jj(iip)=ipk; vv(iip)=(Tn(yc,iky)*Tn(xc,ikx)); iip=iip+1;

    end
end
% g(ip)=g(ip)-1;

ix=nx; iy=1; ip=iy+(ix-1)*ny;
for ikx=0:nt-1
    for iky=0:nt-1
        ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
        g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
%         J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
        ii(iip)=ip; jj(iip)=ipk; vv(iip)=(Tn(yc,iky)*Tn(xc,ikx)); iip=iip+1;

    end
end
g(ip)=g(ip)-1;

ix=nx; iy=ny; ip=iy+(ix-1)*ny;
for ikx=0:nt-1
    for iky=0:nt-1
        ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
        g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
%         J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
        ii(iip)=ip; jj(iip)=ipk; vv(iip)=(Tn(yc,iky)*Tn(xc,ikx)); iip=iip+1;

    end
end
g(ip)=g(ip)-1;

J=J+sparse([ii(1:iip-1);nx*ny],[jj(1:iip-1);nx*ny],[vv(1:iip-1);0]);

% toc;
% tic;
%
%bords
ii=zeros(nt*nt*nx*ny,1); jj=ii; vv=ii; iip=1;
%left

for iy=1+1:length(y)-1
    ix=1; ip=iy+(ix-1)*ny;
    for ikx=0:nt-1
        for iky=0:nt-1
            ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
            g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
%             J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
ii(iip)=ip; jj(iip)=ipk; vv(iip)=(Tn(yc,iky)*Tn(xc,ikx)); iip=iip+1;

        end
    end
%     g(ip)=g(ip)-1; % bc gives 2nd order maybe
%     g(ip)=g(ip)-(y(iy)>0); % bc gives 2nd order maybe
%     g(ip)=g(ip)-(1-y(iy).^2); % bc maybe some 12th order, really smooth, in the center point
%     g(ip)=g(ip)-(1-abs(y(iy))); % bc maybe 3rd
%     g(ip)=g(ip)-exp(-1./(1-y(iy).^2)); % bc

    %right
    ix=nx; ip=iy+(ix-1)*ny;
    for ikx=0:nt-1
        for iky=0:nt-1
            ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
            g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
%             J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
            ii(iip)=ip; jj(iip)=ipk; vv(iip)=(Tn(yc,iky)*Tn(xc,ikx)); iip=iip+1;

        end
    end
      g(ip)=g(ip)-1; %
%       g(ip)=g(ip)-(1-y(iy))/2; %
end

%bottom
for ix=1+1:length(x)-1
    iy=1; ip=iy+(ix-1)*ny;
    for ikx=0:nt-1
        for iky=0:nt-1
            ipk=iky+1+(ikx-1+1)*ny; xc=x(ix); yc=y(iy);
            g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
%             J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
ii(iip)=ip; jj(iip)=ipk; vv(iip)=(Tn(yc,iky)*Tn(xc,ikx)); iip=iip+1;

        end
    end

%     g(ip)=g(ip)-1; % bc
%       g(ip)=g(ip)-(1+x(ix))/2; %
      g(ip)=g(ip)-exp((-1+x(ix))/2/0.03); %

%top
    iy=ny; ip=iy+(ix-1)*ny;
    for ikx=0:nt-1
        for iky=0:nt-1
            ipk=iky+1+(ikx-1+1)*ny;xc=x(ix); yc=y(iy);
            g(ip)=g(ip)+u(ipk)*(Tn(yc,iky)*Tn(xc,ikx));
%             J(ip,ipk)=J(ip,ipk)+(Tn(yc,iky)*Tn(xc,ikx));
ii(iip)=ip; jj(iip)=ipk; vv(iip)=(Tn(yc,iky)*Tn(xc,ikx)); iip=iip+1;

        end
    end
      g(ip)=g(ip)-(1+x(ix))/2; %
end
J=J+sparse([ii(1:iip-1);nx*ny],[jj(1:iip-1);nx*ny],[vv(1:iip-1);0]);

% toc;
%
% tic;
u=u-J\g;
% toc;
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
                uecc=uecc+u(ipk)*(Tn(0.1,iky)*Tn(0.1,ikx));
            end
        end


%
ix=(length(x)+1)/2; iy=(length(y)+1)/2; ip=iy+(ix-1)*ny; 
% uarr=[uarr;uPhys(ip)]
uarr=[uarr;uecc];
narr=[narr; nt];
[narr,uarr]
end

%%

% clf; 
[ev,evs]=eig(full(J),full(B)); evs=diag(evs); 
% [ev,evs]=eigs((J),(B),20,"smallestabs"); evs=diag(evs); 
% plot(evs,'x');
% aa=load("evsFD.mat"); hold on;
% plot(aa.evs,'o');
% uarr=[uarr;real(evs(1:5)')]


narr=[narr; nt];

%%
save('evsCheb',"uarr","narr");
%%
mesh(reshape(log10(abs(u)),[ny,nx]))
% mesh(reshape((abs(u)),[ny,nx]))
%%

clf;
% semilogy(narr,abs(uarr-2/pi/pi),'x-'); hold on;
% semilogy(narr,abs(uarr-uarr(end)),'x-'); hold on;
loglog(narr,abs(uarr-uarr(end,:)),'x-'); hold on;
% semilogy(narr,narr.^(-2),'o-')
loglog(narr,narr.^(-2),'o-')
loglog(narr,narr.^(-3),'o-')
loglog(narr,narr.^(-4),'o-')
loglog(narr,narr.^(-5),'o-')
loglog(narr,narr.^(-12),'o-')

%%
d=diff(uarr)
d(1)/d(2)
%%
clf; 
[ev,evs]=eig(full(J),full(B)); evs=diag(evs); 
% [ev,evs]=eigs((J),(B),20,"smallestabs"); evs=diag(evs); 
plot(evs,'x');
% aa=load("evsFD.mat"); hold on;
% plot(aa.evs,'o');

%hold on; text(real(evs),imag(evs),num2str([1:length(evs)]'));
%%
iev=1;
%%
% x=-cos(linspace(0,pi,nx*5))'; y=x;
x=-cos(linspace(0,pi,nx))'; y=x;
% u=real(ev(:,iev)); iev=iev+1;
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
% title("ev:"+num2str(evs(iev)))
%%
uarr=[uarr;uPhys((length(x)+1)/2)]
%%
d=diff(uarr)
d(1)/d(2)
%%
close all;
plot(x,uPhys,'-x')

%%
clf;
plot(dTn2(x,n)); hold on; 
% plot(-n^2*cos(n*acos(x))./sin(n*acos(x)).^2+n*sin(n*acos(x)).*cos(acos(x))./sin(acos(x)).^3);
plot(n*((n+1)*cos(n*acos(x))-sin((n+1)*acos(x))./sin(acos(x)))./(x.^2-1))

%%
clf;
loglog(narr,abs(uarr-2/pi/pi),'x-'); hold on;
loglog(narr,narr.^(-2),'-')
loglog(narr,narr.^(-4),'-')

%%
save('chebcenter','narr','uarr');
