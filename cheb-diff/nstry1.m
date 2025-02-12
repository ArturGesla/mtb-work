clc; clear;
N=3; %GL points, polynomial of order at most N-1
R=10; H=1; 

%Grids - fully staggered (Canuto p. 152)
% xu=-cos(linspace(0,pi,N));
ru=[-cos(pi*(0:N-1)/(N-1)),1]; %GL points
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

rp=rv; 
zp=zv;

%%
u=zeoros((N+1)^2*4); % as always p u v w 
J=sparse((N+1)^2*4,(N+1)^2*4,0);
g=u; % conti umom vmom wmom

%conti du/dr *2/R + u *2/(R(r+1)) + dw/dz *2/H =0
for ix=2:N
    for iy=2:N
        xc=rp(ix); zc=zp(iy);
ip=((iy-1)+(ix-1)*(N+1))*4+1;

    end
end



