clc; clf;
x=-10:0.1:10; y=-10:0.1:10; uG=0; vG=0;
[X,Y]=meshgrid(x,y); arsc=3; alw=4;


alpha=0*pi/180; xc=0; yc=0; 
Xvec=[reshape(X,length(x)*length(y),1),reshape(Y,length(x)*length(y),1)];
Xvec=Xvec*[cos(alpha), -sin(alpha); sin(alpha),cos(alpha)];
Xr=reshape(Xvec(:,1),length(x),length(y));
Yr=reshape(Xvec(:,2),length(x),length(y));
phi=1/2*(Xr.^2-Yr.^2);
[u,v]=gradient(phi,mean(diff(x)),mean(diff(x)));
r=sqrt(X.^2+Y.^2); u=u.*exp(-r.^2/1); v=v.*exp(-r.^2/1); 
X1=X+xc; Y1=Y+yc; uG=uG+interp2(X1,Y1,u,X,Y,"linear",0); vG=vG+interp2(X1,Y1,v,X,Y,"linear",0);
quiver(xc,yc,[cos(alpha)],[sin(alpha)],arsc/2,Alignment="tail",Color='r',LineWidth=alw,MaxHeadSize=0.8); hold on;
quiver(xc,yc,-[cos(alpha)],[sin(alpha)],arsc/2,Alignment="tail",Color='r',LineWidth=alw,MaxHeadSize=0.8); 
quiver(xc,yc,[-sin(alpha)],[cos(alpha)],arsc/2,Alignment="head",Color='r',LineWidth=alw,MaxHeadSize=0.8); 
quiver(xc,yc,[-sin(alpha)],-[cos(alpha)],arsc/2,Alignment="head",Color='r',LineWidth=alw,MaxHeadSize=0.8); 
%

alpha=-65*pi/180; xc=4; yc=0.7; 
Xvec=[reshape(X,length(x)*length(y),1),reshape(Y,length(x)*length(y),1)];
Xvec=Xvec*[cos(alpha), -sin(alpha); sin(alpha),cos(alpha)];
Xr=reshape(Xvec(:,1),length(x),length(y));
Yr=reshape(Xvec(:,2),length(x),length(y));
phi=1/2*(Xr.^2-Yr.^2);
[u,v]=gradient(phi,mean(diff(x)),mean(diff(x)));
r=sqrt(X.^2+Y.^2); u=u.*exp(-r.^2/1); v=v.*exp(-r.^2/1); 
X1=X+xc; Y1=Y+yc; uG=uG+interp2(X1,Y1,u,X,Y,"linear",0); vG=vG+interp2(X1,Y1,v,X,Y,"linear",0);
quiver(xc,yc,[cos(alpha)],[sin(alpha)],arsc/2,Alignment="tail",Color='r',LineWidth=alw,MaxHeadSize=0.8); hold on;
quiver(xc,yc,-[cos(alpha)],-[sin(alpha)],arsc/2,Alignment="tail",Color='r',LineWidth=alw,MaxHeadSize=0.8); 
quiver(xc,yc,[-sin(alpha)],[cos(alpha)],arsc/2,Alignment="head",Color='r',LineWidth=alw,MaxHeadSize=0.8); 
quiver(xc,yc,-[-sin(alpha)],-[cos(alpha)],arsc/2,Alignment="head",Color='r',LineWidth=alw,MaxHeadSize=0.8); 

%
alpha=35*pi/180; xc=1.5; yc=5.2; 
Xvec=[reshape(X,length(x)*length(y),1),reshape(Y,length(x)*length(y),1)];
Xvec=Xvec*[cos(alpha), -sin(alpha); sin(alpha),cos(alpha)];
Xr=reshape(Xvec(:,1),length(x),length(y));
Yr=reshape(Xvec(:,2),length(x),length(y));
phi=1/2*(Xr.^2-Yr.^2);
[u,v]=gradient(phi,mean(diff(x)),mean(diff(x)));
r=sqrt(X.^2+Y.^2); u=u.*exp(-r.^2/1); v=v.*exp(-r.^2/1); 
X1=X+xc; Y1=Y+yc; uG=uG+interp2(X1,Y1,u,X,Y,"linear",0); vG=vG+interp2(X1,Y1,v,X,Y,"linear",0);
quiver(xc,yc,[cos(alpha)],[sin(alpha)],arsc/2,Alignment="tail",Color='r',LineWidth=alw,MaxHeadSize=0.8); hold on;
quiver(xc,yc,-[cos(alpha)],-[sin(alpha)],arsc/2,Alignment="tail",Color='r',LineWidth=alw,MaxHeadSize=0.8); 
quiver(xc,yc,[-sin(alpha)],[cos(alpha)],arsc/2,Alignment="head",Color='r',LineWidth=alw,MaxHeadSize=0.8); 
quiver(xc,yc,-[-sin(alpha)],-[cos(alpha)],arsc/2,Alignment="head",Color='r',LineWidth=alw,MaxHeadSize=0.8); 

%

alpha=100*pi/180; xc=-1.5; yc=3.3; 
Xvec=[reshape(X,length(x)*length(y),1),reshape(Y,length(x)*length(y),1)];
Xvec=Xvec*[cos(alpha), -sin(alpha); sin(alpha),cos(alpha)];
Xr=reshape(Xvec(:,1),length(x),length(y));
Yr=reshape(Xvec(:,2),length(x),length(y));
phi=1/2*(Xr.^2-Yr.^2);
[u,v]=gradient(phi,mean(diff(x)),mean(diff(x)));
r=sqrt(X.^2+Y.^2); u=u.*exp(-r.^2/1); v=v.*exp(-r.^2/1); 
X1=X+xc; Y1=Y+yc; uG=uG+interp2(X1,Y1,u,X,Y,"linear",0); vG=vG+interp2(X1,Y1,v,X,Y,"linear",0);
quiver(xc,yc,[cos(alpha)],[sin(alpha)],arsc/2,Alignment="tail",Color='r',LineWidth=alw,MaxHeadSize=0.8); hold on;
quiver(xc,yc,-[cos(alpha)],-[sin(alpha)],arsc/2,Alignment="tail",Color='r',LineWidth=alw,MaxHeadSize=0.8); 
quiver(xc,yc,[-sin(alpha)],[cos(alpha)],arsc/2,Alignment="head",Color='r',LineWidth=alw,MaxHeadSize=0.8); 
quiver(xc,yc,-[-sin(alpha)],-[cos(alpha)],arsc/2,Alignment="head",Color='r',LineWidth=alw,MaxHeadSize=0.8); 

%

% mesh(X,Y,uG)
% quiver(X,Y,uG,vG); axis equal;
% %%
% psi=cumtrapz(Y(:,1),u,1)-cumtrapz(X(1,:),v,2);
% contour(X,Y,psi)
% % quiver(X,Y,u,v)

%


f= @(t,x) [interp2(X,Y,uG,x(1),x(2)); interp2(X,Y,vG,x(1),x(2))];
[t,y]=rk4_2(@(t,x)f(t,x),2500,0.1,0,[0.1;-2]);
%
 % quiver(X,Y,uG,vG);
hold on; plot(y(1,:),y(2,:),'k-',LineWidth=alw); axis equal;
quiver(y(1,end),y(2,end),diff(y(1,end-1:end)),diff(y(2,end-1:end)),arsc*1000,Color='k',LineWidth=alw,MaxHeadSize=10);

mar=1.5; 
minx=min(y(1,:));  maxx=max(y(1,:)); xlim([(maxx+minx)/2-(maxx-minx)/2*mar,(maxx+minx)/2+(maxx-minx)/2*mar])
minx=min(y(2,:));  maxx=max(y(2,:)); ylim([(maxx+minx)/2-(maxx-minx)/2*mar,(maxx+minx)/2+(maxx-minx)/2*mar])
ylim([-2.5 7])
xlim([-3 5.5])

text(-2.4,0.5,"saddle point","Interpreter","latex",'Color','red',"FontSize",22)
text(0.5,-1,"turbulent trajectory","Interpreter","latex",'Color','black',"FontSize",22,Rotation=0)

fnts=10; jfm_plt_aid_comm;
axis off;
%% 
exportgraphics(gcf,"turbulenceSketch.eps")

%% Hopf
clf;
xd=-0.6:0.1:0; xu=0:0.01:0.9;

plot(xd,xd*0,'k',LineWidth=1); hold on;
plot(xu,xu*0,'k--',LineWidth=1);
plot(xu,sqrt(xu)*0.9,'k-',LineWidth=1);
text(0,-0.1,"$\mu_c$",'Interpreter','latex',FontSize=15)
text(-0.85,0.9,"$A$",'Interpreter','latex',FontSize=15)
plot([-0.6 -0.6],[0 1],'k',LineWidth=1); quiver(-0.6,0.9,0,0.15,Color='black',MaxHeadSize=20,LineWidth=1)

axis off; 
axis equal;
fnts=10; jfm_plt_aid_comm; size_sq23;
set(gcf,"Position",[   450   361   247   130])

exportgraphics(gcf,"hopf1_meca.eps")

%% h1
clf;
mu=0.2;
f= @(t,x) [mu*x(1)-x(2)-x(3)*x(1);mu*x(2)+x(1)-x(3)*x(2);-x(3)+x(1)*x(1)+x(2)*x(2)];
[t,y]=rk4_2(@(t,x)f(t,x),400,0.1,0,[0;0.1/2;0]);
% plot(y')

plot(y(1,:),y(2,:),'k'); hold on; axis equal; set(gca,"XTick",[0]);set(gca,"YTick",[0]); set(gca,"XTickLabels",[]);set(gca,"YTickLabels",[]); axis square;
xlabel("x"); ylabel("y"); title("$\mu>\mu_c$"); grid on;

fnts=10; jfm_plt_aid_comm; size_sq23;
set(gcf,"Position",[     450   337   163   154])
exportgraphics(gcf,"hopf2_meca.eps")

clf;
mu=-0.09;
f= @(t,x) [mu*x(1)-x(2)-x(3)*x(1);mu*x(2)+x(1)-x(3)*x(2);-x(3)+x(1)*x(1)+x(2)*x(2)];
[t,y]=rk4_2(@(t,x)f(t,x),400,0.1,0,[0;0.1/2;0]);
% plot(y')

plot(y(1,:),y(2,:),'k'); hold on; axis equal; set(gca,"XTick",[0]);set(gca,"YTick",[0]); set(gca,"XTickLabels",[]);set(gca,"YTickLabels",[]); axis square;
xlabel("x"); ylabel("y"); title("$\mu<\mu_c$"); grid on;

fnts=10; jfm_plt_aid_comm; size_sq23;
set(gcf,"Position",[     450   337   163   154])
exportgraphics(gcf,"hopf3_meca.eps")

