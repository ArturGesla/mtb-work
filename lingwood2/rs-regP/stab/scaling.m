clc; clear;
cd C:\Users\Artur\Documents\rotst\plots_yd

Re=500;
a=importdata("u"+num2str(Re)+"-600-160.dat");
% a=importdata("u"+num2str(Re)+".dat.base");
a=loadrs2(a,600,160);
 cd C:\Users\Artur\Documents\GitHub\mtb-work\lingwood2\rs-regP\stab;
%%


omiar=[];
xa=[0.15;    0.0];
% xa=[0.25;0.06]; irng=200:400;
xa=[0.25;0.01]; irng=150:350; %Reh=500



for    iii=    irng


R=a.xc(1,iii)*sqrt(Re);
bbar=0;

zc=a.zc(:,1);
zw=(zc(1:end-1)+zc(2:end))/2;
zw(end+1)=2*zw(end)-zw(end-1);

 x=xa(:,end);
xa=[x];
eps=1e-6;

%profile mapping

data.x=(zw*sqrt(Re))';
xc=a.xc(1,iii);
ir=iii;
u=mean(a.u(:,ir-1:ir)')'; u=u/xc;
p=(a.p(:,ir)); p(2:end-1)=(p(2:end-1)-p(2))*Re;
v=(a.v(:,ir)); v=v/xc;
w=(a.w(:,ir)); w=w*sqrt(Re);

uu=reshape([p,u,v,w]',[162*4,1]);
data.u=uu;



zh=@(x,y) imag(imagOmega(x,y,bbar,R,data));

%%

for ii=1:10
d2fdx2=(zh(x(1)+eps,x(2))-2*zh(x(1),x(2))+zh(x(1)-eps,x(2)))/eps^2;
d2fdxdy=(zh(x(1)+eps,x(2)+eps)+zh(x(1)-eps,x(2)-eps)-zh(x(1)+eps,x(2)-eps)-zh(x(1)-eps,x(2)+eps))/eps^2/4;
d2fdy2=(zh(x(1),x(2)+eps)-2*zh(x(1),x(2))+zh(x(1),x(2)-eps))/eps^2;

jac=[d2fdx2,d2fdxdy; d2fdxdy, d2fdy2];
g=[(zh(x(1)+eps,x(2))-zh(x(1)-eps,x(2)))/2/eps;(zh(x(1),x(2)+eps)-zh(x(1),x(2)-eps))/2/eps];
x=x-jac\g;
xa=[xa,x];
det(jac);
norm(g);
% imagOmega(x(1),x(2),bbar,R,data)

fprintf("ar: %4.2e\tai: %4.2e\tR: %4.2f\tnorm(g): %4.2e\tomi: %4.2e\n",x(1),x(2),R,norm(g),zh(x(1),x(2)));
if(norm(g)<1e-9), break; end
end
if(norm(g)>1e-7), error("no conv"); end

omiar(end+1)=zh(x(1),x(2));
end

Rar=a.xc(1,irng)*sqrt(Re);
%%

save('omiar2D500.mat','Rar','omiar')
%%

a=load('omiar2D500.mat'); plot(a.Rar,a.omiar,'-x'); grid on; hold on;
a=load('omiarsup50Reh500.mat'); plot(a.Rar,a.omiar,'r-'); grid on; hold on;
% a=load('omiarbelow40.mat'); plot(a.Rar,a.omiar,'r-'); grid on; hold on;

xlabel("R"); ylabel("abs growth rate");
legend("2d profile substitution","1d","1d")
title("om0i self s vs 2d Re="+num2str(Re))

%%


ir=100;

u=mean(a.u(:,ir-1:ir)')';
p=(a.p(:,ir));
v=(a.v(:,ir));
w=(a.w(:,ir));
zc=a.zc(:,ir);
xc=a.xc(1,ir);
r=xc*sqrt(Re);

clf;
plot(u/xc,zc*sqrt(Re),'-o'); hold on; grid on;
plot(v/xc,zc*sqrt(Re),'-o'); 
plot(w*sqrt(Re),zc*sqrt(Re),'-o'); 
plot((p(2:end-1)-p(2))*Re,zc(2:end-1)*sqrt(Re),'-o'); 

cd C:\Users\Artur\Documents\GitHub\mtb-work\lingwood2\rs-regP\stab;
a2=load("../rs-np-162-k-0.313-reh-1000.mat");
up=reshape(a2.u,[4,length(a2.u)/4])';
x=a2.x;
z=x;
zw=[2*z(1)-z(2),z(1:end-1),2*z(end-1)-z(end-2)];
zc=(zw(1:end-1)+zw(2:end))/2;
zw=zw(2:end);
plot(up(:,1:end),zc,'k+-'); 




