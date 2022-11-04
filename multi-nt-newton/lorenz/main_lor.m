clc; clear; close all;
%
neq=3; r=28; b=8/3; sigma=10;

T=1.5586;
%
np=400; %valid points, no repeats
dt=T/(np);
g=zeros(neq*np+1,1);
J=zeros(neq*np+1);
u=zeros(neq*np+1,1)+T;

main_lorenz_ti

x=X(:,1); x=x+rand(1,length(x))'.*x*0.1;
y=X(:,2); y=y+rand(1,length(y))'.*y*0.1;
z=X(:,3); y=y+rand(1,length(y))'.*y*0.1;
u(1:neq:end-1)=x(1:end-1); u(2:neq:end-1)=y(1:end-1); u(3:neq:end-1)=z(1:end-1);
% plot(x,y); axis equal;

uM=[];
uM=[uM,u];
%%
% calc J and g
 for i=1:8
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
    g(iz+1)=(u(ix)-u(ixpp))/ds/2-60.9946;
    J(iz+1,ix)=1/ds/2;
    J(iz+1,ixpp)=-1/ds/2;

    %
    du=-J\g;
   
    fprintf('%s\n',"iter: "+num2str(i)+" res: "+num2str(norm(du)))
    u=u+du;
    uM=[uM,u];
 end
%%
hold on;
% u=u-du2;
plot(u(1:neq:end-1),u(2:neq:end-1))
% plot(u(1:2:end)+utang(1:2:end),u(2:2:end)+utang(2:2:end))
% plot([u(1:2:end), u(1:2:end)+utang(1:2:end)]',[u(2:2:end),u(2:2:end)+utang(2:2:end)]'); axis equal
plot([u(1:neq:end-1), u(1:neq:end-1)+du(1:neq:end-1)]',[u(2:neq:end-1),u(2:neq:end-1)+du(2:neq:end-1)]','-o'); axis equal
% plot([u(1:2:end), u(1:2:end)+du(1:2:end)]',[u(2:2:end),u(2:2:end)+du(2:2:end)]','-x'); 
grid on; grid minor;

%%
hold on;
plot3(uM(1:neq:end-1,:),uM(2:neq:end-1,:),uM(3:neq:end-1,:))
% plot3([uM(1:neq:end-1,end-1), uM(1:neq:end-1,end-1)+du(1:neq:end-1)]',[uM(2:neq:end-1,end-1),uM(2:neq:end-1,end-1)+du(2:neq:end-1)]',[uM(3:neq:end-1,end-1), uM(3:neq:end-1,end-1)+du(3:neq:end-1)]','-o'); 
%%
spy(J); grid on; grid minor;
%%
plot(uM(1:neq:end-1,:))
%%
close all;
plot(uM(1:neq:end-1,:),uM(2:neq:end-1,:))
%%
close all;
f=figure(Position=[2200 202 911 598]); fnts=14;
% f=figure(Position=[2200 202 300 300]); fnts=14;
% hold on;
set(f,'defaulttextinterpreter','latex')


hold on;
plot(uM(1:2:end-1,1),uM(2:2:end-1,1))
plot(uM(1:2:end-1,end),uM(2:2:end-1,end),'-o')

x=r*cos(linspace(0,2*pi,np+1));
y=r*sin(linspace(0,2*pi,np+1));
u(1:2:end-1)=x(1:end-1); u(2:2:end-1)=y(1:end-1);
plot(u(1:2:end-1),u(2:2:end-1),'-x')
axis equal;
% grid on;

h1=legend("iteration=0","iteration=5, resid=1e-16","analytical solution",Location="best");
set(h1, 'Interpreter','latex')
grid on; grid minor; 
set(gca,"FontSize",fnts,"FontName","Latin Modern Math");
xlabel("x"); ylabel("y"); title("Multi timestep newton method $ \vert$ subHoph $\mu=-0.015 $ $ \vert$ 100 points for period")
exportgraphics(gcf,'plot.png','Resolution',200)

%%
hold on;
plot(uM(1,:)',uM(2,:)'); axis equal;
plot(u(1:2:end),u(2:2:end))
%%
plot(sqrt(u(1:2:end).^2+u(2:2:end).^2))
%%
hold on;
plot(uM')
plot(sqrt(sum(uM.^2,1))')
