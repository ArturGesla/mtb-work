% Noack system
% f=[mu*u-v*g-w*u;
%     mu*v+u/g-v*w;
%     -w+u^2+v^2*g*g];

clc; clear; close all;
%
%parms
% mu=0.04; r=sqrt(mu); gm=1; %gamma
mu=8; r=sqrt(mu); gm=1; %gamma
neq=3; np=100; np=np+2;

%init
t=0:2*pi/(np-1):2*pi; u=r*cos(t); v=r*sin(t)/gm; w=r^2;
T=2*pi; x=u; y=v; z=w;

% np=np-1;
dt=T/(np-1);
g=zeros(neq*np+1,1);
J=sparse(neq*np+1,neq*np+1);%zeros(neq*np+1);
u=zeros(neq*np+1,1)+T;

u(1:neq:end-1)=x; u(2:neq:end-1)=y; u(3:neq:end-1)=z;
%
% %
% % plot(x,y); axis equal;
% [a,bb]=(min(abs(gradient(x(1:end-1))))); 
% %  bb=np-1; 
% phaseIndex=bb;%phase index
% ds=1/np; derX=0;%(x(b-1)-x(b+1))/2/ds
% 
% uM=[];
% uM=[uM,u];
% uMC=[];
% r=15;
%
% calc J and g
% r=r*1.1
% for ii=1:42
    for ii=1:1
tic;
 for i=1:15
    
     evalJacRhs
%      evalJacRhs_bdf2
    %
    du=-sparse(J)\g;
    
%    [L,U] = ilu(sparse(J),struct('type','ilutp','droptol',1e-16));
% du=bicgstab(-J,g,1e-10,6,L,U);
% du=bicgstab(-J,g,1e-10,60);
   
    fprintf('%s\n',"iter: "+num2str(i)+" res: "+num2str(norm(du)))
    u=u+du;
%     uM=[uM,u];
 end
 toc;
%  uMC=[uMC,u];
%  r=r+(24.74-r)*0.1
% r=r-0.1
    end
%
close all;
X=reshape(u(1:end-4),[3,np-1])'; 
plot(X(:,1),X(:,2)); hold on; plot(X(1,1),X(1,2),'o'); plot(X(2,1),X(2,2),'>'); axis equal; grid on;

% stab
Jc=full(J(1:end-1,1:end-1)); B=zeros(length(Jc)); B(end-2:end,1:3)=eye(3); Jc(end-2:end,1:3)=0;
[evc,evs]=eig(Jc,B); evs=diag(evs); b=evs<Inf; evs=evs(b);
%
u=evc(:,1)*0
u=evc(:,1)*0; u(end-2:end)=1e-3;
v=Jc\u
x3=reshape(v(1:end-3),[3,np-1])';
%
close all;

x1=reshape(evc(1:end-3,1),[3,np-1])'; 
x2=reshape(evc(1:end-3,2),[3,np-1])'; 
% x3=reshape(evc(1:end-3,3),[3,np-1])'; 

xp=X; plot3(xp(:,1),xp(:,2),xp(:,3)); hold on; grid on;
xp=X+x1; plot3(xp(:,1),xp(:,2),xp(:,3)); hold on; grid on;
xp=X+x2+x1; plot3(xp(:,1),xp(:,2),xp(:,3)); hold on; grid on;
xp=X+x3; plot3(xp(:,1),xp(:,2),xp(:,3),'-x'); hold on; grid on;
% view([0,0,1]); 
legend(["base";num2str(evs)])

%% f analysis
Z=fft(X(:,1));
plot(abs(Z)); hold on; 


z1=fft(x1(:,1));
plot(abs(z1),'x-'); hold on; 

grid on; %set(gca,"yscale","log");
set(gca,"xscale","log"); 