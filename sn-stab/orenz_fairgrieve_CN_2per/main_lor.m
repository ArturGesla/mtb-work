clc; clear; close all;
%
% neq=3; r=28; b=8/3; sigma=10; T=1.5586; np=70; %valid points, no repeats
neq=3; r=24; b=8/3; sigma=10; np=101; np=np+2;%valid points, no repeats
% neq=3; r=160; b=8/3; sigma=10; T=1.1521; np=240; %valid points, no repeats
% neq=3; r=140; b=8/3; sigma=10; T=1.5586;
% neq=3; r=145; b=8/3; sigma=10; T=1.5586;
% r=15;
%
main_lorenz_ti; 
%
close all;
plot(X(:,1),X(:,2),'x-'); hold on;
plot(X(1,1),X(1,2),'o-')

%
% np=np;
% dt=T/(np-1);
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

uM=[];
uM=[uM,u];
uMC=[];
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
    uM=[uM,u];
 end
 toc;
 uMC=[uMC,u];
%  r=r+(24.74-r)*0.1
% r=r-0.1
end
%%
hold on; close all;
% u=u-du2;
plot(u(1:neq:end-1),u(2:neq:end-1),'-x'); hold on;
plot(u(end-1-2),u(end-1-1),'-o')
plot(u(1),u(2),'-sq')
% plot(u(1:2:end)+utang(1:2:end),u(2:2:end)+utang(2:2:end))
% plot([u(1:2:end), u(1:2:end)+utang(1:2:end)]',[u(2:2:end),u(2:2:end)+utang(2:2:end)]'); axis equal
% plot([u(1:neq:end-1), u(1:neq:end-1)+du(1:neq:end-1)]',[u(2:neq:end-1),u(2:neq:end-1)+du(2:neq:end-1)]','-o'); axis equal
% plot([u(1:2:end), u(1:2:end)+du(1:2:end)]',[u(2:2:end),u(2:2:end)+du(2:2:end)]','-x'); 
grid on; grid minor; axis equal;
%% cum visu
% plot(uMC(1:neq:end-1,:),uMC(2:neq:end-1,:));
plot(uM(1:neq:end-1,:),uM(2:neq:end-1,:));
xlabel("x"); ylabel("y"); title("Lorenz | UPO near Hopf at r=24.73")
%% Fairgrieve part
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
lam=1./(1-evs) % one of lambda should be 1
exp=1./T/2*log(abs(lam))
% %%
% v=evc(:,2);
% x=reshape(v,[3,np])';
% X=reshape(u(1:end-1),[3,np])';
% 
% %% deco
% C=A(1:end-3,1:end-3); D=A(1:end-3,end-2:end);
% E=A(end-2:end,1:end-3); F=A(end-2:end,end-2:end);
% G=(-E*inv(C)*D+F);
% eigs()

%% visu 
close all;
hold on;
% u=u-du2;
axis equal;
plot(u(1:neq:end-1),u(2:neq:end-1)); 
% plot(u(1:neq:3),u(2:neq:2*neq),'o');
plot(u(end-3),u(end-2),'sq');
%
% ie=2;
% plot(u(1:2:end-1)+evc(1:2:end-1,ie),u(2:2:end-1)+evc(2:2:end-1,ie));
% evc=evc*20;
ie=3;
mult=10;
plot(u(1:neq:end-1)+real(evc(1:neq:end,ie)*mult),u(2:neq:end-1)+real(evc(2:neq:end,ie))*mult);
plot([u(1:neq:end-1),u(1:neq:end-1)+real(evc(1:neq:end,ie)*mult)]',[u(2:neq:end-1), u(2:neq:end-1)+real(evc(2:neq:end,ie)*mult)]');
title("Lorenz | r: "+num2str(r)+" fl mult: "+num2str(lam(ie)))

