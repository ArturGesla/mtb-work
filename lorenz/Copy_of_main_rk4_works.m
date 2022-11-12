% clc; clear;
sigma=10;
b=8/3;
r=24;

c=sqrt(b*(r-1));
x=c;
y=c;
z=r-1;
en=x.^2+y^2+z^2
% x=0;
% y=0;
% z=0;

J=[-sigma, sigma,0;r-z,-1,-x;y,x,-b];

[v,u]=eig(J);
u=diag(u);

u
v
evc=real(v(:,2));
% alpha=7;
% a0=6; a1=7;
% alpha=(a0+a1)/2;
% alpha=   5.906250000000000
% alpha=   6.144226074218750 % dt=0.05
alpha=   5.894348144531250 % dt=0.05/2/2


% 
% alpha=5
a=alpha
x=x+alpha*evc(1);
y=y+alpha*evc(2);
z=z+alpha*evc(3);

%
% r=28;
% sigma=10;
% b=8/3;
nt=8e3;
dt=0.05/2/2;

[T,X]=rk4(@odefun,nt,dt,0,[x;y;z],r,b,sigma);

% [Tt,Y]=ode45(@(t,y)odefun(t,y,r,b,sigma),[0 nt*dt], [x,y,z]);

hold on;
plot(T,X); T=T'; X=X';
% plot(Tt,Y,'o');
%%
close all;
plot3(X(:,1),X(:,2),X(:,3),'-')
%% analysis
z=fft(X(:,1));
f=0:1/T(end):1/mean(diff(T));
close all; plot(f(2:end),abs(z(2:end))); grid on; grid minor; 
[a,d]=max(abs(z(2:end)));
per=1/f(d+1)
close all; semilogy(f(2:end)./f(d+1),abs(z(2:end))); grid on; grid minor; 
np=per/dt
%%
x=X(end,1);
y=X(end,2);
z=X(end,3);
%%
[T,X]=rk4(@odefun,round(np)+1,dt,0,[x;y;z],r,b,sigma);
plot(T,X); T=T'; X=X';
% plot(Tt,Y,'o');
%%
close all;
% plot3(X(:,1),X(:,2),X(:,3),'-x')
plot(X(:,1),X(:,2),'-x')
x=X(1:end-1,1); x=x';
y=X(1:end-1,2); y=y';
z=X(1:end-1,3); z=z';

%% interp
% t=linspace(0,per,np*2+1);
% x=interp1(T,X(:,1),t);
% y=interp1(T,X(:,2),t);
% z=interp1(T,X(:,3),t);
% 
% close all;hold on;
% plot3(x,y,z,'o'); 
% plot3(X(:,1),X(:,2),X(:,3),'-x')
% 
% x=x(1:end-1);
% y=y(1:end-1);
% z=z(1:end-1);
%% mtsnm
neq=3; %r=28; b=8/3; sigma=10;
% r=15;
% T=1.5586;
T=per;
%
np=size(x,2); %valid points, no repeats
dt=T/(np);
g=zeros(neq*np+1,1);
J=sparse(neq*np+1,neq*np+1);%zeros(neq*np+1);
u=zeros(neq*np+1,1)+T;

% main_lorenz_ti

% x=X(:,1); %x=x+rand(1,length(x))'.*x*0.1;
% y=X(:,2); %y=y+rand(1,length(y))'.*y*0.1;
% z=X(:,3); %y=y+rand(1,length(y))'.*z*0.1;
u(1:neq:end-1)=x(1:end); u(2:neq:end-1)=y(1:end); u(3:neq:end-1)=z(1:end);
% plot(x,y); axis equal;

uM=[];
uM=[uM,u];
uMC=[];
res=[];
rM=[];
% r=15;
%%
% plot3(u(1:neq:end-1),u(2:neq:end-1),u(3:neq:end-1),'x-')
close all;
plot(u(1:neq:end-1),u(2:neq:end-1),'x-')
%%
for ic=1:90
r=r-0.1
tic;
 for i=1:15
    
     for ip=1:np
        ix=ip*neq-2;
        iy=ip*neq-1;
        iz=ip*neq;

        ixp=mod(ip*neq-2+neq-1,neq*np)+1;
        ixpp=mod(ip*neq-2+neq+neq-1,neq*np)+1;
        iyp=mod(ip*neq-1+neq-1,neq*np)+1;
        izp=mod(ip*neq-0+neq-1,neq*np)+1;
        ixm=mod(ip*neq-2-neq-1,neq*np)+1;
        iym=mod(ip*neq-1-neq-1,neq*np)+1;
        izm=mod(ip*neq-0-neq-1,neq*np)+1;

%     x=u(ix); y=u(iy); z=u(iz);% r=sqrt(x^2+y^2);
% 
%     xNext=u(ixp); yNext=u(iyp); zNext=u(izp);
%     xPrev=u(ixm); yPrev=u(iym); zPrev=u(izm);

%     xNextNext=u(mod(ip*2-1+neq+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);

    T=u(neq*np+1); dt=T/np; ds=1/np;

    g(ix)=T*(sigma*(u(iy)-u(ix)))           +(u(ixm)-u(ixp))/ds/2;
    g(iy)=T*(r*u(ix)-u(iy)-u(ix)*u(iz))     +(u(iym)-u(iyp))/ds/2;
    g(iz)=T*(u(ix)*u(iy)-b*u(iz))           +(u(izm)-u(izp))/ds/2;

    
    
    J(ix,ix)=T*(sigma*(-1));
    J(ix,iy)=T*(sigma*(1));
    J(ix,ixp)=(-1)/ds/2;
    J(ix,ixm)=(1)/ds/2;
    J(ix,neq*np+1)=(sigma*(u(iy)-u(ix)));

    J(iy,iy)=T*(-1);
    J(iy,ix)=T*(r-u(iz));
    J(iy,iz)=T*(-u(ix));
    J(iy,iyp)=(-1)/ds/2;
    J(iy,iym)=(1)/ds/2;
    J(iy,neq*np+1)=(r*u(ix)-u(iy)-u(ix)*u(iz));

    J(iz,iz)=T*(-b);
    J(iz,ix)=T*(u(iy));
    J(iz,iy)=T*(u(ix));
    J(iz,izp)=(-1)/ds/2;
    J(iz,izm)=(1)/ds/2;
    J(iz,neq*np+1)=(u(ix)*u(iy)-b*u(iz));
  
    end

    % last row
%     g(iz+1)=(u(ix)-u(ixpp))/ds/2-60.9946;
%     J(iz+1,ix)=1/ds/2;
%     J(iz+1,ixpp)=-1/ds/2;
    
    g(iz+1)=(u(ixm)-u(ixp))/ds/2+13.212372967119705;
    J(iz+1,ixm)=1/ds/2;
    J(iz+1,ixp)=-1/ds/2;

    %
    du=-sparse(J)\g;
    
%    [L,U] = ilu(sparse(J),struct('type','ilutp','droptol',1e-16));
% du=bicgstab(-J,g,1e-10,6,L,U);
% du=bicgstab(-J,g,1e-10,60);
   
    fprintf('%s\n',"iter: "+num2str(i)+" res: "+num2str(norm(du)))
    u=u+du;
    uM=[uM,u];
 end
 toc;
 uMC=[uMC,u];
 res=[res;norm(du)];
 rM=[rM,r];
end
 %%
 plot3(uM(1:neq:end-1,:),uM(2:neq:end-1,:),uM(3:neq:end-1,:))
 %%
 plot3(uMC(1:neq:end-1,:),uMC(2:neq:end-1,:),uMC(3:neq:end-1,:))
 %%
 plot(uMC(1:neq:end-1,:),uMC(2:neq:end-1,:))
 %%
 close all;
 plot3(uM(1:neq:end-1,[1,end]),uM(2:neq:end-1,[1,end]),uM(3:neq:end-1,[1,end]))