clc; clear; close all;
%
% neq=3; r=28; b=8/3; sigma=10; T=1.5586; np=70; %valid points, no repeats
neq=3; r=24; b=8/3; sigma=10; np=1000; %valid points, no repeats
% neq=3; r=160; b=8/3; sigma=10; T=1.1521; np=240; %valid points, no repeats
% neq=3; r=140; b=8/3; sigma=10; T=1.5586;
% neq=3; r=145; b=8/3; sigma=10; T=1.5586;
% r=15;
%
main_lorenz_ti

dt=T/(np);
g=zeros(neq*np+1,1);
J=sparse(neq*np+1,neq*np+1);%zeros(neq*np+1);
u=zeros(neq*np+1,1)+T;

gstab=zeros(neq*np,1);
Jstab=sparse(neq*np,neq*np);%zeros(neq*np+1);
ustab=zeros(neq*np,1);


amp=0;
x=X(:,1); x=x+rand(1,length(x))'.*x*amp;
y=X(:,2); y=y+rand(1,length(y))'.*y*amp;
z=X(:,3); y=y+rand(1,length(y))'.*y*amp;
u(1:neq:end-1)=x(1:end-1); u(2:neq:end-1)=y(1:end-1); u(3:neq:end-1)=z(1:end-1);
% plot(x,y); axis equal;
[a,bb]=(min(abs(gradient(x(1:end-1))))); 
%  bb=np-1; 
phaseIndex=bb;%phase index
ds=1/np; derX=0;%(x(b-1)-x(b+1))/2/ds

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
% plot([u(1:neq:end-1), u(1:neq:end-1)+du(1:neq:end-1)]',[u(2:neq:end-1),u(2:neq:end-1)+du(2:neq:end-1)]','-o'); axis equal
% plot([u(1:2:end), u(1:2:end)+du(1:2:end)]',[u(2:2:end),u(2:2:end)+du(2:2:end)]','-x'); 
grid on; grid minor; axis equal;
%
n=reshape(reshape(gstab,[neq,np])./vecnorm(reshape(gstab,[neq,np])),[neq*np,1]);
P=eye(length(gstab))-n*n';
mask=[];
for i=1:np
mask=blkdiag(mask,ones(neq,neq));
end
P=P.*mask;
%

evM=[];
for i=1:np
A=P*Jstab*P;
% A=P*Jstab;
ip=i;  A=A(1+(ip-1)*neq:3+(ip-1)*neq,1+(ip-1)*neq:3+(ip-1)*neq);
[evc,ev]=eig(A); ev=diag(ev); %plot(real(ev),imag(ev)+i,'o'); grid on; hold on;
evM=[evM;ev'];
end

for i=1:np
ip=i;index=1+(ip-1)*neq:3+(ip-1)*neq;
if(max(evM(i,:)>1e-6))
plot(u(index(1)),u(index(2)),'rx');
else
plot(u(index(1)),u(index(2)),'bx');
end
end
% xlim([-80 80]); ylim([-80 80]);
xlabel("x"), ylabel("y"); title("Lorenz r="+num2str(r)+"; b=8/3; sigma=10");
exportgraphics(gcf,"lorenzstabR"+num2str(r)+".png","resolution",150)
%% stab
% n=gstab/norm(gstab);
n=reshape(reshape(gstab,[neq,np])./vecnorm(reshape(gstab,[neq,np])),[neq*np,1]);
P=eye(length(gstab))-n*n';
mask=[];
for i=1:np
mask=blkdiag(mask,ones(neq,neq));
end
P=P.*mask;
%
close all;
A=P*Jstab*P;
ip=5; ip=ip+1;
% A=A(1+ip*neq:3+ip*neq,1+ip*neq:3+ip*neq);
% A=Jstab;
% A=1/2*(Jstab+transpose(Jstab));
% A=J;
% A=J(1:end-1,1:end-1);
% B=eye(length(u));
% B(end,end)=0;
% ev=eig(A,B); plot(ev,'x'); grid on; hold on;
[evc,ev]=eig(A); ev=diag(ev); plot(real(ev),imag(ev),'o'); grid on; hold on;
%%
iev=42;
du=evc(:,iev)*3;
% plot([u(1:neq:end-1), u(1:neq:end-1)+du(1:neq:end-1)]',[u(2:neq:end-1),u(2:neq:end-1)+du(2:neq:end-1)]','-o'); axis equal
plot3([u(1:neq:end-1), u(1:neq:end-1)+du(1:neq:end-1)]',[u(2:neq:end-1),u(2:neq:end-1)+du(2:neq:end-1)]',[u(3:neq:end-1),u(3:neq:end-1)+du(3:neq:end)]','-o'); 
axis equal; grid on;
title(num2str(ev(iev)))
%%
evM=[];
for i=1:np
A=P*Jstab*P;
% A=P*Jstab;
ip=i;  A=A(1+(ip-1)*neq:3+(ip-1)*neq,1+(ip-1)*neq:3+(ip-1)*neq);
[evc,ev]=eig(A); ev=diag(ev); plot(real(ev),imag(ev)+i,'o'); grid on; hold on;
evM=[evM;ev'];
end
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

%% maybe global
clc; close all;
% JJ=eye(np*neq)+dt*A;
JJ=P+dt*A;

Jglob=eye(neq);
for i=1:np
    Jglob=(JJ(1+(i-1)*neq:neq+(i-1)*neq,1+(i-1)*neq:neq+(i-1)*neq))*Jglob;
end
eig(Jglob)
% exp(u(end)*max(ev))