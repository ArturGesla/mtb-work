clc; clear;
N=11; %GL points, polynomial of order at most N-1
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
u(3:4:end)=av;
u(2:4:end)=au;
u(1:4:end)=ap;
u(4:4:end)=aw;
dlmwrite("u.dat",u,'precision',17);
system(" cp u.dat /home/gesla/Documents/git/rotst2/build");


save('input');
%%

tic;
% [g,jac]=evalRhsAndJac(rp,zp,ru,zu,rv,zv,rw,zw,N,R,H,u); 
evalRhsAndJac(); 
% !time ./run_evalRhsAndJac.sh /home/gesla/mtb-installation-folder/
toc;

%
a=load("output.mat");
g=a.g;
%
% mesh(reshape(g(1:4:end),[N+1,N+1])) % conti
mesh(reshape(g(1:4:end),[N+1,N+1]))
%%
mesh(reshape(log10(abs(au)),[N+1,N+1]))


%%
cd /home/gesla/Documents/git/rotst2/build
a=importdata("g.dat"); n=12;
mesh(reshape(a(2:4:end),[n,n]))