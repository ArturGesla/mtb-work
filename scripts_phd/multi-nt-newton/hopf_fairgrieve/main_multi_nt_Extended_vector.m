clc; clear; close all;
%
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
np=50;
dt=T/(np-1);
g=zeros(2*np+1,1);
gstab=zeros(2*np,1);
J=zeros(2*np+1);
Jstab=zeros(2*np);
% J2=zeros(2*np);
u=zeros(2*np+1,1)+T;
% r=sqrt(-mu)
% r=0.05

phi0=0;%pi/8;
x=r*cos(linspace(phi0+0,phi0+2*pi,np)); %x=x+rand(1,length(x)).*x*0.1;
y=r*sin(linspace(phi0+0,phi0+2*pi,np)); %y=y+rand(1,length(y)).*y*0.1;
%
% x=r*cos(linspace(0,2*pi,np+1)); %x=x+rand(1,length(x)).*x*0.1;
% y=r*sin(linspace(0,2*pi,np+1)); %y=y+rand(1,length(y)).*y*0.1;
%
u(1:2:end-1)=x(1:end); u(2:2:end-1)=y(1:end);
% plot(x,y); axis equal;

uM=[];
uM=[uM,u];
%
% calc J and g
 for i=1:10
    for ip=1:np-1
    x=u(ip*2-1); y=u(ip*2); r=sqrt(x^2+y^2);
    ixp=mod(ip*2-1+neq-1,neq*np)+1; iyp=ixp+1;
    ixm=mod(ip*2-1-neq-1,neq*np)+1-neq*((ip*2-1-neq-1)<0); iym=ixm+1;
    
    xNext=u(ixp); yNext=u(iyp);
    xPrev=u(ixm); yPrev=u(iym);
    
    T=u(2*np+1); %dt=T/np; 
    ds=1/(np-1);
%     T=T; dt=T/np; ds=1/np;
    g(ip*2-1)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2))+(xPrev-xNext)/ds/2;
    g(ip*2)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2))+(yPrev-yNext)/ds/2;
    
    
    
    J(ip*2-1,ip*2-1)=T*((mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x);
    J(ip*2-1,ip*2)=T*(x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y));
    J(ip*2-1,ixp)=-1/ds/2;
    J(ip*2-1,ixm)=1/ds/2;
  
    J(ip*2,ip*2-1)=T*(y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x));
    J(ip*2,ip*2)=T*((mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y);
    J(ip*2,iyp)=-1/ds/2;
    J(ip*2,iym)=1/ds/2;

    %last column
    J(ip*2-1,end)=((mu+r^2-r^4)*x-y*(omega+b*r^2));
    J(ip*2,end)=((mu+r^2-r^4)*y+x*(omega+b*r^2));
  
    %jstab
% T=1;
% gstab(ip*2-1)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2));
%     gstab(ip*2)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2));
% 
%     Jstab(ip*2-1,ip*2-1)=T*((mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x);
%     Jstab(ip*2-1,ip*2)=T*(x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y));
% % %     J(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/ds/2;
% %     J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=1/ds/2;
%   
%     Jstab(ip*2,ip*2-1)=T*(y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x));
%     Jstab(ip*2,ip*2)=T*((mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y);
% %     J(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/ds/2;
% %     J(ip*2,mod(ip*2-neq-1,neq*np)+1)=1/ds/2;
% 
% %     %last column
% %     J(ip*2-1,end)=((mu+r^2-r^4)*x-y*(omega+b*r^2));
% %     J(ip*2,end)=((mu+r^2-r^4)*y+x*(omega+b*r^2));
% %   
%     
    end

    ip=np;
    x=u(ip*2-1); y=u(ip*2); r=sqrt(x^2+y^2);
    ix=ip*2-1; iy=ix+1;
    ixp=mod(ip*2-1+neq-1,neq*np)+1; iyp=ixp+1;
    ixm=mod(ip*2-1-neq-1,neq*np)+1-neq*((ip*2-1-neq-1)<0); iym=ixm+1;
    
%     xNext=u(mod(ip*2-1+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);
%     xNextNext=u(mod(ip*2-1+neq+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);
%     xPrev=u(mod(ip*2-1-neq-1,neq*np)+1); yPrev=u(mod(ip*2-neq-1,neq*np)+1);
    T=u(2*np+1); 
%     dt=T/np; 
ds=1/(np-1);
%     T=T; dt=T/np; ds=1/np;
%     g(ix)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2))+(u(ixm)-u(ixp))/ds/2;
%     g(iy)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2))+(u(iym)-u(iyp))/ds/2;

    g(ix)=u(ix)-u(1);
    g(iy)=u(iy)-u(2);
    
  J(ix,ix)=1; J(ix,1)=-1; 
  J(iy,iy)=1; J(iy,2)=-1; 
    
%     J(ip*2-1,ip*2-1)=T*((mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x);
%     J(ip*2-1,ip*2)=T*(x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y));
%     J(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/ds/2;
%     J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=1/ds/2;
%   
%     J(ip*2,ip*2-1)=T*(y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x));
%     J(ip*2,ip*2)=T*((mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y);
%     J(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/ds/2;
%     J(ip*2,mod(ip*2-neq-1,neq*np)+1)=1/ds/2;

    %last column
%     J(ip*2-1,2*np+1)=((mu+r^2-r^4)*x-y*(omega+b*r^2));
%     J(ip*2,2*np+1)=((mu+r^2-r^4)*y+x*(omega+b*r^2));
    
%     %jstab
% T=1;
%   gstab(ip*2-1)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2));
%     gstab(ip*2)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2));
%     Jstab(ip*2-1,ip*2-1)=T*((mu+r^2-r^4)+x*(2*x-2*r^2*2*x)-y*b*2*x);
%     Jstab(ip*2-1,ip*2)=T*(x*(2*y-2*r^2*2*y)-(omega+b*r^2)-y*(2*b*y));
% % %     J(ip*2-1,mod(ip*2-1+neq-1,neq*np)+1)=-1/ds/2;
% %     J(ip*2-1,mod(ip*2-1-neq-1,neq*np)+1)=1/ds/2;
%   
%     Jstab(ip*2,ip*2-1)=T*(y*(2*x-2*r^2*2*x)+(omega+b*r^2)+x*(2*b*x));
%     Jstab(ip*2,ip*2)=T*((mu+r^2-r^4)+y*(2*y-2*r^2*2*y)+x*b*2*y);
% %     J(ip*2,mod(ip*2+neq-1,neq*np)+1)=-1/ds/2;
% %     J(ip*2,mod(ip*2-neq-1,neq*np)+1)=1/ds/2;
% 

    % last row
%     g(2*np+1)=(x-xNextNext)/ds/2;
%     J(2*np+1,mod(ip*2-1+neq*2-1,neq*np)+1)=-1/ds/2;
%     J(2*np+1,mod(ip*2-1-1,neq*np)+1)=1/ds/2;
    
    g(2*np+1)=(u(ix)-u(ixm))/ds;
    J(2*np+1,ix)=1/ds;
    J(2*np+1,ixm)=-1/ds;
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
%%
hold on;
% u=u-du2;
axis equal;
plot(u(1:2:end-1),u(2:2:end-1)); 
plot(u(1),u(2),'o');
plot(u(end-2),u(end-1),'sq');
% u=u1; plot(u(1:2:end-1),u(2:2:end-1))
% u=u2; plot(u(1:2:end-1),u(2:2:end-1))

% plot(u(1:2:end)+utang(1:2:end),u(2:2:end)+utang(2:2:end))
% plot([u(1:2:end), u(1:2:end)+utang(1:2:end)]',[u(2:2:end),u(2:2:end)+utang(2:2:end)]'); axis equal
% plot([u(1:2:end), u(1:2:end)+du(1:2:end)]',[u(2:2:end),u(2:2:end)+du(2:2:end)]','-o'); axis equal
% plot([u(1:2:end), u(1:2:end)+du(1:2:end)]',[u(2:2:end),u(2:2:end)+du(2:2:end)]','-x'); 
grid on; grid minor; axis equal;
xlabel("x"), ylabel("y"); title("subHopf mu=-0.05; b=-1; omega=1"); %legend("lambda = -1.6944","lambda = 0.0944")
% exportgraphics(gcf,"subHopfstab11.png","resolution",150);
% close all;

%% new stab - first ti
 
% g(ip*2-1)=T*((mu+r^2-r^4)*x-y*(omega+b*r^2))+(xPrev-xNext)/ds/2;
%     g(ip*2)=T*((mu+r^2-r^4)*y+x*(omega+b*r^2))+(yPrev-yNext)/ds/2;
%     
% subHopf = @(t,x) [(mu+(x(1)^2+x(2)^2)^2-(x(1)^2+x(2)^2)^4)*x(1)-x(2)*(omega+b*(x(1)^2+x(2)^2)^2);
%     (mu+(x(1)^2+x(2)^2)^2-(x(1)^2+x(2)^2)^4)*x(2)+x(1)*(omega+b*(x(1)^2+x(2)^2)^2)]; % Anonymous Function
 [t,X] = rk4(@subhopfp, 10*np,T/(np-1),mu,b,om,0,[r;0],[0;0]);
% T=1/0.868; [t,X] = ode45(lorenz, [0 :T/(np):1*T], [-5.6775    2.1838  134.4854]);
% plot3(X(:,1),X(:,2),X(:,3))
X=X';
close all;
hold on; 
plot(X(:,1),X(:,2))
plot(X(1,1),X(1,2),'-o')
plot(X(end,1),X(end,2),'-sq')
axis equal; axis square;

%%
close all
plot(t,sqrt(X(:,1).^2+X(:,2).^2)); hold on;
plot(t,(t+1)./(t+1)*r)

%% floquet
close all;
rr=sqrt(X(:,1).^2+X(:,2).^2);
semilogy(t,abs(rr-r)); hold on;
%%
close all;
rl=log(abs(rr-r)); b=abs(rl)<inf; b(1:1000)=0;
a=polyfit(t(b),rl(b),1);
plot(t,rl); hold on;
plot(t,a(1)*t+a(2))

%% Fairgrieve part
A=J;
B(np*neq,np*neq)=1; B(np*neq-1,np*neq-1)=1;  B(np*neq+1,np*neq+1)=0; B=sparse(B);
evs=eigs(A,B,neq,"smallestabs")
% lam=1./(1-evs) % one of lambda should be 1
% exp=1./T*log(abs(lam))


%% old stab
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