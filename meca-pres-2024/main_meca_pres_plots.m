addpath     'C:\Users\Izabela\Documents\git\rotst2\scripts\source_for_mtb'
cd     'C:\Users\Izabela\Documents\git\mtb-work\meca-pres-2024'

%% hopf ti

nP=1; r=24; np=200;
lorenz = @(t,x) [10*(x(2)-x(1));r*x(1)-1*x(1)*x(3)-x(2);x(1)*x(2)-(8/3)*x(3)]; % Anonymous Function
% [T,X] = ode45(lorenz, [0 :T/(np):1*T], [-6.2262, -11.0027,13.0515]);
ev=[ -0.0328   -0.6228    0.6818]*5;
% ev=[ -0.0328   -0.6228    0.6818]*-2;
T=0.6779; T= 0.6803*1.2; [t,X] = ode45(lorenz, [0 :T/(np-2):T/(np-2)*(np-1)*nP],[10.135982315094342  10.189521543725682  25.691556187487929]+ev);
T=0.6779; T= 0.6803*1.5; [t,X2] = ode45(lorenz, [0 :T/(np-2):T/(np-2)*(np-1)*nP],[10.135982315094342  10.189521543725682  25.691556187487929]+ev*0);
% T=0.6779; T= 0.6803; [t,X] = ode45(lorenz, [0 :T/(np-2):T/(np-2)*(np-1)*nP], [9.3 11.1 22.2]);
%
close all;
plot(X(:,1),X(:,2)); hold on;
plot(X2(:,1),X2(:,2));
xlabel("x"); ylabel("y"); grid on;

fnts=10; jfm_plt_aid_comm; size_sq23;
%%
% exportgraphics(gcf,"hopforb24.eps")
exportgraphics(gcf,"hopforb24-shoot.eps")

%%
a=load("udat.mat"); np=(length(a.u)-1)/3;
X=reshape(a.u(1:end-1),[3,np])';

close all;
plot(X(:,1),X(:,2)); hold on;
plot(X(:,1),X(:,2),'x'); hold on;
% plot(X2(:,1),X2(:,2));
xlabel("x"); ylabel("y"); grid on;

fnts=10; jfm_plt_aid_comm; size_sq23;
%%
exportgraphics(gcf,"hopforb24-coll.eps")

%% lorenz sn

close all;

% up=[evc(:,105);0]; ntp=nt; u1=up(1:end-1);
% up=[j2*u2;0]; ntp=nt; 
u0=load("uSN-lorenz-nt15.mat"); u0=u0.u;
up=[u0]; ntp=nt; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
zp=[zp(1:2,:);zeros(1000,3)];
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp); xpp=xp;
% xp=xp+xpb;  
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xBase=xp;
plot(xp(:,1),xp(:,2)); grid on; hold on; xBase=xp;
xlabel("x"); ylabel("y"); grid on;

plot(X2(:,1),X2(:,2)); legend("nt = 1","exact","Location","best")

fnts=10; jfm_plt_aid_comm; size_sq23;
%%
exportgraphics(gcf,"hopforb24-sn-1.eps")
%% harmonics
close all;


% up=[evc(:,105);0]; ntp=nt; u1=up(1:end-1);
% up=[j2*u2;0]; ntp=nt; 
u0=load("uSN-lorenz-nt15.mat"); u0=u0.u;
up=[u0]; ntp=nt; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
zp=[zp(1:1,:)*0;zp(2,:);zeros(1000,3)];
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp); xpp=xp;
% xp=xp+xpb;  
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xBase=xp;
plot(linspace(0,1,length(xp)),xp(:,1)); grid on; hold on; xBase=xp;
% xlabel("t/T"); 
ylabel("x"); grid on;

% plot(X2(:,1),X2(:,2)); legend("nt = 1","exact","Location","best")

fnts=10; jfm_plt_aid_comm; size_sq23; set(gcf,"Position",[450.0000  553.0000  431.0000   88.3333]);
% set(gca,"Position",[0.1472    0.3521    0.7578    0.5729])
%%
exportgraphics(gcf,"hopforb24-sn-har-1.eps")

%% langford
l=1.99; close all;
mu=l;

delta=0.8*mu-0.8*2.8+1;
z=[-1-sqrt(delta),-1+sqrt(delta)]/-0.4;
r=sqrt(-z.*(z-mu));

langford = @(t,y) [(l-3)*y(1)-1/4*y(2)+y(1)*(y(3)+0.2*(1-y(3)^2));
    (l-3)*y(2)+1/4*y(1)+y(2)*(y(3)+0.2*(1-y(3)^2));
    l*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))]; % Anonymous Function
 % [t,X] = ode45(langford, [0,1000], [ 0.2334 -0.9401 1.0858]);
 [t,X0] = ode45(langford, [0,1000], [ r(2),0,z(2)]);
% plot(X(:,1),X(:,2))
% plot3(X(:,1),X(:,2),X(:,3)); hold on;

l=2.005;
langford = @(t,y) [(l-3)*y(1)-1/4*y(2)+y(1)*(y(3)+0.2*(1-y(3)^2));
    (l-3)*y(2)+1/4*y(1)+y(2)*(y(3)+0.2*(1-y(3)^2));
    l*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))]; % Anonymous Function
 % [t,X] = ode45(langford, [0,1000], [ 0.2334 -0.9401 1.0858]);
 [t,X1] = ode45(langford, [0,1000], [0.2493   -1.0089    1.3466]);

 plot3(X0(:,1),X0(:,2),X0(:,3),LineWidth=2); hold on;
 plot3(X1(:,1),X1(:,2),X1(:,3)); hold on;

 grid; xlabel('x'); ylabel('y'); zlabel('z'); legend("$\lambda=1.99$","$\lambda=2.005$","Location","northwest")
 fnts=10; jfm_plt_aid_comm; size_sq23;

 %%
 exportgraphics(gcf,"langford1.png","Resolution",300)

 %% lang cheb
 load lang.mat;

 close all;
up=u+1e-3*sum([evc(:,1:2);0,0],2)*0;
neq=3;
xch=X*0;
for i=0:nt-1    
    xch(:,1)=xch(:,1)+up(i*neq+1).*cos(i*acos(tch1));
    xch(:,2)=xch(:,2)+up(i*neq+2).*cos(i*acos(tch1));
    xch(:,3)=xch(:,3)+up(i*neq+3).*cos(i*acos(tch1));
end

% plot(tch,xch,'-'); hold on; set(gca,"ColorOrderIndex",1); %same as ycut
plot3(xch(:,1),xch(:,2),xch(:,3),'-'); hold on; %cd mecaplot3(xch(1,1),xch(1,2),xch(1,3),'>'); hold on; 
grid on; hold on;

up=u+1e-2*sum([evc(:,1:2);0,0],2);
neq=3;
xch=X*0;
for i=0:nt-1    
    xch(:,1)=xch(:,1)+up(i*neq+1).*cos(i*acos(tch1));
    xch(:,2)=xch(:,2)+up(i*neq+2).*cos(i*acos(tch1));
    xch(:,3)=xch(:,3)+up(i*neq+3).*cos(i*acos(tch1));
end

% plot(tch,xch,'-'); hold on; set(gca,"ColorOrderIndex",1); %same as ycut
plot3(xch(:,1),xch(:,2),xch(:,3),'-'); hold on; %plot3(xch(1,1),xch(1,2),xch(1,3),'>'); hold on; 
grid on; hold on;

grid; xlabel('x'); ylabel('y'); zlabel('z'); %legend("$\lambda=1.99$","$\lambda=2.005$","Location","northwest")
 fnts=10; jfm_plt_aid_comm; size_sq23; grid;
%%
exportgraphics(gcf,"langford2.png","Resolution",300)