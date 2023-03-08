clc; clear; close all;
%%
neq=2; mu=-0.05; b=-1; omega=1;

% sub
% mu=-0.1;
p = [-1 0 1 0 mu 0];
r = roots(p);
r2=r(3);
r1=r(5);
r=r1;
om=1;
% b=1;
omEff=om+b*r^2;
f=omEff/2/pi;
T=1/f;
%
np=250;
dt=T/np;
g=zeros(2*np+1,1);
gstab=zeros(2*np,1);
J=zeros(2*np+1);
Jstab=zeros(2*np);
% J2=zeros(2*np);
u=zeros(2*np+1,1)+T;
% r=sqrt(-mu)
% r=0.05

phi0=0;%pi/8;
x=r*cos(linspace(phi0+0,phi0+2*pi,np+1)); %x=x+rand(1,length(x)).*x*0.1;
y=r*sin(linspace(phi0+0,phi0+2*pi,np+1)); %y=y+rand(1,length(y)).*y*0.1;
% 
% x=r*cos(linspace(0,2*pi,np+1)); %x=x+rand(1,length(x)).*x*0.1;
% y=r*sin(linspace(0,2*pi,np+1)); %y=y+rand(1,length(y)).*y*0.1;
%
u(1:2:end-1)=x(1:end-1); u(2:2:end-1)=y(1:end-1);
% plot(x,y); axis equal;

uM=[];
uM=[uM,u];
%
% calc J and g
 for i=1:10
    for ip=1:np-1
    x=u(ip*2-1); y=u(ip*2); r=sqrt(x^2+y^2);
    xNext=u(mod(ip*2-1+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);
    xPrev=u(mod(ip*2-1-neq-1,neq*np)+1); yPrev=u(mod(ip*2-neq-1,neq*np)+1);
    T=u(2*np+1); dt=T/np; ds=1/np;
%     T=T; dt=T/np; ds=1/np;
    g(ip*2-1)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2))+(xPrev-xNext)/ds/2;
    g(ip*2)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2))+(yPrev-yNext)/ds/2;
    
    
    
    J(ip*2-1,ip*2-1)=T*((mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x);
    J(ip*2-1,ip*2)=T*(x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y));
    J(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/ds/2;
    J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=1/ds/2;
  
    J(ip*2,ip*2-1)=T*(y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x));
    J(ip*2,ip*2)=T*((mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y);
    J(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/ds/2;
    J(ip*2,mod(ip*2-neq-1,neq*np)+1)=1/ds/2;

    %last column
    J(ip*2-1,end)=((mu+r^2-r^4)*x-y*(omega+b*r^2));
    J(ip*2,end)=((mu+r^2-r^4)*y+x*(omega+b*r^2));
  
    %jstab
T=1;
gstab(ip*2-1)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2));
    gstab(ip*2)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2));

    Jstab(ip*2-1,ip*2-1)=T*((mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x);
    Jstab(ip*2-1,ip*2)=T*(x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y));
% %     J(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/ds/2;
%     J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=1/ds/2;
  
    Jstab(ip*2,ip*2-1)=T*(y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x));
    Jstab(ip*2,ip*2)=T*((mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y);
%     J(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/ds/2;
%     J(ip*2,mod(ip*2-neq-1,neq*np)+1)=1/ds/2;

%     %last column
%     J(ip*2-1,end)=((mu+r^2-r^4)*x-y*(omega+b*r^2));
%     J(ip*2,end)=((mu+r^2-r^4)*y+x*(omega+b*r^2));
%   
    
    end

    ip=np;
    x=u(ip*2-1); y=u(ip*2); r=sqrt(x^2+y^2);
    xNext=u(mod(ip*2-1+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);
    xNextNext=u(mod(ip*2-1+neq+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);
    xPrev=u(mod(ip*2-1-neq-1,neq*np)+1); yPrev=u(mod(ip*2-neq-1,neq*np)+1);
    T=u(2*np+1); dt=T/np; ds=1/np;
%     T=T; dt=T/np; ds=1/np;
    g(ip*2-1)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2))+(xPrev-xNext)/ds/2;
    g(ip*2)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2))+(yPrev-yNext)/ds/2;
    
   
    
    J(ip*2-1,ip*2-1)=T*((mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x);
    J(ip*2-1,ip*2)=T*(x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y));
    J(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/ds/2;
    J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=1/ds/2;
  
    J(ip*2,ip*2-1)=T*(y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x));
    J(ip*2,ip*2)=T*((mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y);
    J(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/ds/2;
    J(ip*2,mod(ip*2-neq-1,neq*np)+1)=1/ds/2;

    %last column
    J(ip*2-1,2*np+1)=((mu+r^2-r^4)*x-y*(omega+b*r^2));
    J(ip*2,2*np+1)=((mu+r^2-r^4)*y+x*(omega+b*r^2));
    
    %jstab
T=1;
  gstab(ip*2-1)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2));
    gstab(ip*2)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2));
    Jstab(ip*2-1,ip*2-1)=T*((mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x);
    Jstab(ip*2-1,ip*2)=T*(x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y));
% %     J(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/ds/2;
%     J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=1/ds/2;
  
    Jstab(ip*2,ip*2-1)=T*(y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x));
    Jstab(ip*2,ip*2)=T*((mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y);
%     J(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/ds/2;
%     J(ip*2,mod(ip*2-neq-1,neq*np)+1)=1/ds/2;


    % last row
    g(2*np+1)=(x-xNextNext)/ds/2;
    J(2*np+1,mod(ip*2-1+neq*2-1,neq*np)+1)=-1/ds/2;
    J(2*np+1,mod(ip*2-1-1,neq*np)+1)=1/ds/2;
%     
%  % last row
%     g(2*np+1)=y;
%     J(2*np+1,2*np)=1;
% %     J(2*np+1,mod(ip*2-neq-1,neq*np)+1)=1/ds/2;



    %
    du=-J\g;
   
    norm(du)
    u=u+du;
    uM=[uM,u];
 end
%
hold on;
% u=u-du2;
axis equal;
plot(u(1:2:end-1),u(2:2:end-1))
% u=u1; plot(u(1:2:end-1),u(2:2:end-1))
% u=u2; plot(u(1:2:end-1),u(2:2:end-1))

% plot(u(1:2:end)+utang(1:2:end),u(2:2:end)+utang(2:2:end))
% plot([u(1:2:end), u(1:2:end)+utang(1:2:end)]',[u(2:2:end),u(2:2:end)+utang(2:2:end)]'); axis equal
% plot([u(1:2:end), u(1:2:end)+du(1:2:end)]',[u(2:2:end),u(2:2:end)+du(2:2:end)]','-o'); axis equal
% plot([u(1:2:end), u(1:2:end)+du(1:2:end)]',[u(2:2:end),u(2:2:end)+du(2:2:end)]','-x'); 
grid on; grid minor; axis equal;
xlabel("x"), ylabel("y"); title("subHopf mu=-0.05; b=-1; omega=1"); legend("lambda = -1.6944","lambda = 0.0944")
exportgraphics(gcf,"subHopfstab11.png","resolution",150);
close all;

%
% n=gstab/norm(gstab);
n=reshape(reshape(gstab,[neq,np])./vecnorm(reshape(gstab,[neq,np])),[neq*np,1]);
P=eye(length(gstab))-n*n';
mask=[];
for i=1:np
mask=blkdiag(mask,ones(2,2));
end
P=P.*mask;

close all;
A=P*Jstab*P';
% A=Jstab;
% A=1/2*(Jstab+transpose(Jstab));
% A=J;
% A=J(1:end-1,1:end-1);
% B=eye(length(u));
% B(end,end)=0;
% ev=eig(A,B); plot(ev,'x'); grid on; hold on;
[evc,ev]=eig(A); ev=diag(ev); plot(real(ev),imag(ev),'o'); grid on; hold on; ev

%

%
% iev=1;
% ip=7;
% % iiev=iev+(ip-1)*2;
% iiev=1;
% % utang=evc(:,);
% utang=evc(:,iiev);
% 
% % ustab=u(1:end-1);
% 
% % utang=Jstab*ustab*0.1;
% 
% hold on;
% axis equal;
% plot(u(1:2:end-1),u(2:2:end-1));
% plot(u(1:2:end-1)+utang(1:2:end),u(2:2:end-1)+utang(2:2:end));
% plot([u(1:2:end-1)';u(1:2:end-1)'+utang(1:2:end)'],[u(2:2:end-1)';u(2:2:end-1)'+utang(2:2:end)']);
% % iiev=38; utang=evc(:,iiev); plot([u(1:2:end-1)';u(1:2:end-1)'+utang(1:2:end)'],[u(2:2:end-1)';u(2:2:end-1)'+utang(2:2:end)']);
% title("ip: "+num2str(ip)+" iev: "+num2str(iev)+" lam: "+num2str(ev(iiev)))

% utang=gstab;
% plot([u(1:2:end-1)';u(1:2:end-1)'+utang(1:2:end)'],[u(2:2:end-1)';u(2:2:end-1)'+utang(2:2:end)']);

% maybe global
clc; close all;
% JJ=eye(np*neq)+dt*A;
JJ=P+dt*A;

Jglob=eye(2);
for i=1:np
    Jglob=(JJ(1+(i-1)*neq:2+(i-1)*neq,1+(i-1)*neq:2+(i-1)*neq))*Jglob;
end
eig(Jglob)
exp(u(end)*max(ev))