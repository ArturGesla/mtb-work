clc; clear;
N=31; %GL points, polynomial of order at most N-1
R=1; H=1;

%Grids - fully staggered (Canuto p. 152)
% xu=-cos(linspace(0,pi,N));
ru=[-cos(pi*(0:N-1)/(N-1))]; %GL points
zu=[-1,-cos(pi*((0:N-2)+0.5)/(N-1)),1]; %G points + border
% y2=-cos((2*(1:N)-1)*pi/2/N)
% y2=-cos((2*(1:(N-1))-1)*pi/2/(N-1))

% clf;
% plot(xu,xu,'-x'); hold on;
% plot(yu,yu,'-o')
% plot(y2,y2,'-<')

rw=zu;
zw=ru;

rv=rw;
zv=zu;

rp=rv(2:end-1);
zp=zv(2:end-1);

%
u=zeros((N+1)^2*4,1); % as always p u v w
% J=sparse((N+1)^2*4,(N+1)^2*4,0);
% J=sparse((N+1)^2*4,(N+1)^2*4,0);
g=u; % conti umom vmom wmom

%
nsInterp;
u(2:4:end)=au;
u(4:4:end)=aw;

%

%conti du/dr *2/R + u *2/(R(r+1)) + dw/dz *2/H =0
for ix=1:length(rp)
    for iy=1:length(zp) %grid on conti
        xc=rp(ix); zc=zp(iy);
        ip=((iy)+(ix)*(N+1))*4+1;

        for ikx=0:length(ru)-1 %expansion of velocity
            for iky=0:length(zu)-1
                ikp=(iky+(ikx)*(N+1))*4+2;

                g(ip)=g(ip)+u(ikp)*dTn(xc,ikx)*Tn(zc,iky)*2/R;
                g(ip)=g(ip)+u(ikp)*Tn(xc,ikx)*Tn(zc,iky)*2/(R*(xc+1));
            end
        end

        for ikx=0:length(rw)-1 %expansion of velocity
            for iky=0:length(zw)-1

                ikp=(iky+(ikx)*(N+1))*4+4;

                g(ip)=g(ip)+u(ikp)*Tn(xc,ikx)*dTn(zc,iky)*2/H;

            end
        end

        %

    end
end

%%
mesh(reshape(g(1:4:end),[N+1,N+1]))


