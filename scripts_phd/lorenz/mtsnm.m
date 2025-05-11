clc; clear; close all;
%
eigen_dir
%%
close all;

neq=3; %r=28; b=8/3; sigma=10;
% r=15;
% T=1.5586;
T=per;
%
np=size(X,1)-1; %valid points, no repeats
dt=T/(np);
g=zeros(neq*np+1,1);
J=sparse(neq*np+1,neq*np+1);%zeros(neq*np+1);
u=zeros(neq*np+1,1)+T;

% main_lorenz_ti

x=X(:,1); %x=x+rand(1,length(x))'.*x*0.1;
y=X(:,2); %y=y+rand(1,length(y))'.*y*0.1;
z=X(:,3); %y=y+rand(1,length(y))'.*z*0.1;
u(1:neq:end-1)=x(1:end-1); u(2:neq:end-1)=y(1:end-1); u(3:neq:end-1)=z(1:end-1);
% plot(x,y); axis equal;

uM=[];
uM=[uM,u];
uMC=[];
% r=15;
%%
% calc J and g
% r=r*1.1
tic;
 for i=1:10
    
     evalJacRhs
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
% plot3(uM(1:neq:end-1,:),uM(2:neq:end-1,:),uM(3:neq:end-1,:))
plot3(uM(1:neq:end-1,1),uM(2:neq:end-1,1),uM(3:neq:end-1,1))
plot3(uM(1:neq:end-1,end),uM(2:neq:end-1,end),uM(3:neq:end-1,end))
% plot3([uM(1:neq:end-1,end-1), uM(1:neq:end-1,end-1)+du(1:neq:end-1)]',[uM(2:neq:end-1,end-1),uM(2:neq:end-1,end-1)+du(2:neq:end-1)]',[uM(3:neq:end-1,end-1), uM(3:neq:end-1,end-1)+du(3:neq:end-1)]','-o'); 
%%
hold on;
% plot3(uM(1:neq:end-1,:),uM(2:neq:end-1,:),uM(3:neq:end-1,:))
plot3(uMC(1:neq:end-1,:),uMC(2:neq:end-1,:),uMC(3:neq:end-1,:))
% plot3(uM(1:neq:end-1,end),uM(2:neq:end-1,end),uM(3:neq:end-1,end))
grid on;
%dla r=15 najmniejsza , 92 najwieksza 
%%
spy(J); grid on; grid minor;
%%
plot(uM(1:neq:end-1,:))
%%
close all;
plot(uM(1:neq:end-1,:),uM(2:neq:end-1,:))
%%
close all;
% f=figure(Position=[2200 202 911 598]); fnts=14;
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

% h1=legend("iteration=0","iteration=5, resid=1e-16","analytical solution",Location="best");
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
