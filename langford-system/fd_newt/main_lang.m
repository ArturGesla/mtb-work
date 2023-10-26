close all; clear; clc;
nt=3000; %generalised does not work wtf xd
dt=1e-1;
t0=0;
% a=-3; b=-9.3; c=8; d=-3; e=5.98; close all;
% a=-3; b=-1; c=-b; d=a; e=2/3; close all; x0=[1e-2;0;1e-3]; %per orbit
 mu=1.90; a=mu-3; b=-1/4; c=-b; d=a; e=mu; c2=0.2; close all; x0=[0.1;0;2]; %per orbit
% a=-3; b=-8; c=8; d=-3; e=5.98; close all; x0=[0.1;0.1;0.1]; %per orbit
% langfordG = @(t,y) [a*y(1)+b*y(2)+y(1)*y(3);
%     c*y(1)+d*y(2)+y(2)*y(3);
%     e*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))];
langfordG = @(t,y) [a*y(1)+b*y(2)+y(1)*y(3)+c2*y(1)*(1-y(3)^2);
                    c*y(1)+d*y(2)+y(2)*y(3)+c2*y(2)*(1-y(3)^2);
                    e*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))];% Anonymous Function
% [t,y]=rk4_2(@(t,y)langfordG(t,y),nt,dt,t0,[1;1;1]);
[t,y]=rk4_2(@(t,y)langfordG(t,y),nt,dt,t0,x0);
% [t,y]=rk4_2(@(t,y)langfordG(t,y),nt,dt,t0,[0.003; 0.001; 0.001]);
% plot3(y(1,:),y(2,:),y(3,:));
% xlabel("x"); ylabel("y");  zlabel("z"); 
%
subplot(1,2,1)
plot(t,y'); hold on;
plot(t,sqrt(y(1,:).^2+y(2,:).^2));
subplot(1,2,2)
plot3(y(1,:),y(2,:),y(3,:));
delta=0.8*mu-0.8*2.8+1
z=[-1-sqrt(delta),-1+sqrt(delta)]/-0.4
r=sqrt(-z.*(z-mu))
sqrt(y(1,end).^2+y(2,end).^2)
y(3,end)
%% anal
% lam=1.5;
% lam=(2-sqrt(1.76))/0.4;
lam=mu;
x=0;y=0; z=lam;
delta=0.8*lam-0.8*2.8+1;
z=(1-sqrt(delta))/0.4;
r=sqrt(-z*(z-lam));
t=6;
x=r*cos(t);y=r*sin(t);

J=[lam-3+z+0.2*(1-z^2), -1/4, x-0.2*z*2*x;
    1/4,lam-3+z+0.2*(1-z^2),  y-0.2*z*2*y;
    -2*x,-2*y,lam-2*z];
eig(J)
%%

neq=3;
% r=24; b=8/3; sigma=10;
np=80; np=np+2;%valid points, no repeats
% neq=3; r=160; b=8/3; sigma=10; T=1.1521; np=240; %valid points, no repeats
% neq=3; r=140; b=8/3; sigma=10; T=1.5586;
% neq=3; r=145; b=8/3; sigma=10; T=1.5586;
% r=15;
%
% main_lorenz_ti; 
%
t=linspace(0,2*pi,np);
X=zeros(np,neq);
X(:,1)=r*cos(t); X(:,2)=r*sin(t); X(:,3)=(X(:,1)+1)./(X(:,1)+1)*z;
close all;

plot(X(:,1),X(:,2),'x-'); hold on;
plot(X(1,1),X(1,2),'o-')

%%
T=2*pi/b;
% np=np-1;
dt=T/(np-1);
g=zeros(neq*np+1,1);
J=sparse(neq*np+1,neq*np+1);%zeros(neq*np+1);
u=zeros(neq*np+1,1)+T;

% gstab=zeros(neq*np,1);
% Jstab=sparse(neq*np,neq*np);%zeros(neq*np+1);
% ustab=zeros(neq*np,1);


amp=0;
x=X(:,1); x=x+rand(1,length(x))'.*x*amp;
y=X(:,2); y=y+rand(1,length(y))'.*y*amp;
z=X(:,3); y=y+rand(1,length(y))'.*y*amp;
u(1:neq:end-1)=x(1:end); u(2:neq:end-1)=y(1:end); u(3:neq:end-1)=z(1:end);
% plot(x,y); axis equal;
[a,bb]=(min(abs(gradient(x(1:end-1))))); 
%  bb=np-1; 
phaseIndex=bb;%phase index
ds=1/np; derX=0;%(x(b-1)-x(b+1))/2/ds
%%
uM=[];
uM=[uM,u];
uMC=[];
mu=1.90; a=mu-3; b=-1/4; c=-b; d=a; e=mu; c2=0.2;
% r=15;
%
% calc J and g
% r=r*1.1
% for ii=1:42
    for ii=1:1
tic;
 for i=1:1
    
     evalJacRhs3
%      evalJacRhs_bdf2
    %
    du=-sparse(J)\g;
    
%    [L,U] = ilu(sparse(J),struct('type','ilutp','droptol',1e-16));
% du=bicgstab(-J,g,1e-10,6,L,U);
% du=bicgstab(-J,g,1e-10,60);
   
    fprintf('%s\n',"iter: "+num2str(i)+" res: "+num2str(norm(g)))
    u=u+du;
    uM=[uM,u];
 end
 toc;
 uMC=[uMC,u];
%  r=r+(24.74-r)*0.1
% r=r-0.1
end
%%