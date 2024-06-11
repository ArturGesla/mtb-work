close all; clear; clc;
dt=1e-1;
nt=300/dt; %generalised does not work wtf xd

t0=0; mu=2.005;
% a=-3; b=-9.3; c=8; d=-3; e=5.98; close all;
% a=-3; b=-1; c=-b; d=a; e=2/3; close all; x0=[1e-2;0;1e-3]; %per orbit
 % mu=1.99;
 a=mu-3; b=-1/4; c=-b; d=a; e=mu; c2=0.2; close all; x0=[0.1;0;2]; %per orbit
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
%%
close all;
ra=sqrt(-z.*(z-mu)); ra=ra(2);
r=sqrt(y(1,:).^2+y(2,:).^2);
semilogy(t,abs(r-ra));
r2=abs(r(500:1500)-ra);
t2=t(500:1500);
hold on; semilogy(t2,r2);
a=polyfit(t2,log(r2),1);
semilogy(t2,exp(a(1)*t2))
semilogy(t2,exp(-0.0624*t2))
%% anal
% lam=1.5;
% lam=(2-sqrt(1.76))/0.4;
% mu=1.9;
mu=1.99;
% mu=2.001;
mu=2.005;
lam=mu;
x=0;y=0; z=lam;
delta=0.8*mu-0.8*2.8+1
z=(1-sqrt(delta))/0.4;
r=sqrt(-z*(z-lam));
t=0;
x=r*cos(t);y=r*sin(t);

J=[lam-3+z+0.2*(1-z^2), b, x-0.2*z*2*x;
    -b,lam-3+z+0.2*(1-z^2),  y-0.2*z*2*y;
    -2*x,-2*y,lam-2*z];
evE=eig(J)
expM=[eig(J)];
n=[0;1;0];
msk=[eye(3)-n*n'];
evE=eig(msk'*J*msk)
%%
% a=lam-3+z+0.2*(1-z^2)
% b=x-0.2*2*x*z
% c=-2*x
% d=lam-2*z
% 
% (d+sqrt(d^2+4*c*b))/2
% 
% 
%  r-sqrt(z*(lam-z))
% 1/0.4*(1-sqrt(0.8*lam-1.24))-z
% 
% mu=1/2*(lam-2*z+sqrt((lam-2*z)^2-8*r*(r-0.4*r*z)))

%%

neq=3;
% r=24; b=8/3; sigma=10;ans
np=400; np=np+1;%valid points, no repeats
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
T=2*pi/-b;
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
%
uM=[];
uM=[uM,u];
uMC=[];
% mu=1.90;
a=mu-3; 
% b=-1/4;
c=-b; d=a; e=mu; c2=0.2;

% r=15;
%
% calc J and g
% r=r*1.1
% for ii=1:42
    for ii=1:1
tic;
 for i=1:10
    
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
% fg
A=full(J(1:end-1,1:end-1)); B=[];
B(np*neq,np*neq)=1; B(np*neq-1,np*neq-1)=1; B(np*neq-2,np*neq-2)=1;
% B(3,np*neq-neq)=1; B(2,np*neq-neq-1)=1; B(1,np*neq-neq-2)=1; 
B(np*neq+1,np*neq+1)=0; B=B(1:end-1,1:end-1); %B=sparse(B);
% [evc,evs]=eigs(A,B,neq,"smallestabs"); evs=diag(evs)
% lam=1./(1-evs) % one of lambda should be 1
% exp=1./T*log(abs(lam))

%
[evc,evs]=eig(A,B); evs=diag(evs); 
% 
b2=abs(evs)<Inf; evs=evs(b2); evc=evc(:,b2); 
lam=1./(1-evs); % one of lambda should be 1


exp=1./T*log((lam))
expM=[expM,exp];
%%
clf; axis equal;
plot3(u(1:3:end-1),u(2:3:end-1),u(3:3:end-1));


grid on; hold on;

iev=3; du=evc(:,iev);
du=evc(:,2)+evc(:,1);
% du=evc(:,2)+evc(:,3);
plot3(u(1)+du(1),u(2)+du(2),u(3)+du(3),'sq');
plot3(u(1:3:end-1)+du(1:3:end),u(2:3:end-1)+du(2:3:end),u(3:3:end-1)+du(3:3:end));
% save("exp-fd-nt-"+num2str(np),"exp");

%%
% plot(u(1:3:end-1));
% hold on;
% plot(du(1:3:end)+u(1:3:end-1));
plot(du(1:3:end))