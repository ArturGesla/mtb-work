clc; close all; clear; mua=[];
%
% cd     '/people/gesla/Documents/git/mtb-work/sn-stab/sn-mod-noack';
%

% Lorenz system
 r=24; nt=50; np=4*nt;
%main_lorenz_ti
%

%  sn=load("../lorenz-sn/xforcheb.mat");
% sn=load("../lorenz-sn/xforcheb24.7368.mat");
% sn=load("../lorenz-sn/xforcheb24.73.mat");
sn=load("../lorenz-sn/xforcheb24.mat");
X=sn.xp; t=sn.t; r=sn.r; T=t(end);



 X1=X; 
tch=cos(0:pi/(np*1-1):pi)'; tch1=tch;

%
tch=(tch+1)/2*t(end); X=interp1(t,X,tch); X3=X*0;

ax=sn.u(4:3:sn.nt*3,1);
bx=sn.u(sn.nt*3+1+3:3:sn.nt*3*2,1);
X3(:,1)=cos([1:sn.nt-1].*tch*sn.u(end))*ax*2-sin([1:sn.nt-1].*tch*sn.u(end))*bx*2+sn.u(1);

ay=sn.u(5:3:sn.nt*3,1);
by=sn.u(sn.nt*3+1+4:3:sn.nt*3*2,1);
X3(:,2)=cos([1:sn.nt-1].*tch*sn.u(end))*ay*2-sin([1:sn.nt-1].*tch*sn.u(end))*by*2+sn.u(2);

az=sn.u(6:3:sn.nt*3,1);
bz=sn.u(sn.nt*3+1+5:3:sn.nt*3*2,1);
X3(:,3)=cos([1:sn.nt-1].*tch*sn.u(end))*az*2-sin([1:sn.nt-1].*tch*sn.u(end))*bz*2+sn.u(3);

X=X3;


y=X; v2=[y;flipud(y(2:end-1,:))]; z=real(fft(v2)./length(v2)); a=z; a(2:end,:)=2*a(2:end,:); %z are cheb coeffs
zcut=z*0; zcut([1:nt,end-nt+2:end],:)=z([1:nt,end-nt+2:end],:); ycut=ifft(zcut).*length(zcut); ycut=ycut(1:end/2+1,:);

plot(t,X1,'--'); hold on; set(gca,"ColorOrderIndex",1);
plot(tch,y,':'); hold on; set(gca,"ColorOrderIndex",1);
plot(tch,ycut,'-'); hold on; set(gca,"ColorOrderIndex",1);
plot(tch,X3,'-x'); hold on; set(gca,"ColorOrderIndex",1);
fprintf("cheb init approx accuracy: %4.2e\n",norm(y-ycut,"fro"));
%
u=[reshape(real(a(1:nt,:).'),[3*nt,1])];
u(nt*3+1)=2*pi/T;
% % plot
% neq=3;
% xch=X*0;
% for i=0:nt-1
%     
%     xch(:,1)=xch(:,1)+u(i*neq+1).*cos(i*acos(tch1));
%     xch(:,2)=xch(:,2)+u(i*neq+2).*cos(i*acos(tch1));
%     xch(:,3)=xch(:,3)+u(i*neq+3).*cos(i*acos(tch1));
% end
% 
% plot(tch,xch,'o-'); hold on; set(gca,"ColorOrderIndex",1); %same as ycut
%
% u=u*0; u(end-3)=1;
%
close all; uinit=u;
semilogy(abs(reshape(u(1:end-1),[3,nt])'),'--'); hold on; grid on; grid minor;
%%
for i=1:15
[g,jac]=calculateRhsAndJac(3,nt,u,r);

fprintf("it: %d \t norm(rhs): %4.2e\n",i,norm(g));

if(norm(g)<1e-12)
    break; 
end
u=u-jac\g;
%     u(4:6)'

end
% pause;
semilogy(abs(reshape(u(1:end-1),[3,nt])'),'-x'); hold on; grid on;

%%
%% visu
close all;
up=u;
neq=3;
xch=X*0;
for i=0:nt-1    
    xch(:,1)=xch(:,1)+up(i*neq+1).*cos(i*acos(tch1));
    xch(:,2)=xch(:,2)+up(i*neq+2).*cos(i*acos(tch1));
    xch(:,3)=xch(:,3)+up(i*neq+3).*cos(i*acos(tch1));
end

% plot(tch,xch,'-'); hold on; set(gca,"ColorOrderIndex",1); %same as ycut
plot3(xch(:,1),xch(:,2),xch(:,3),'-'); hold on; plot3(xch(1,1),xch(1,2),xch(1,3),'>'); hold on; 
grid on; hold on;

%

% up=u+[V(:,5);0];
up=uinit;
neq=3;
xch=X*0;
for i=0:nt-1    
    xch(:,1)=xch(:,1)+up(i*neq+1).*cos(i*acos(tch1));
    xch(:,2)=xch(:,2)+up(i*neq+2).*cos(i*acos(tch1));
    xch(:,3)=xch(:,3)+up(i*neq+3).*cos(i*acos(tch1));
end

% plot(tch,xch,'-'); hold on; set(gca,"ColorOrderIndex",1); %same as ycut
plot3(xch(:,1),xch(:,2),xch(:,3),'-'); hold on; plot3(xch(1,1),xch(1,2),xch(1,3),'>'); hold on; 
grid on; hold on;
