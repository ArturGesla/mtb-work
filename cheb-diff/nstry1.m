clc; clear;
N=11; %GL points, polynomial of order at most N-1
R=1; H=1;

%Grids - fully staggered (Canuto p. 152)
% xu=-cos(linspace(0,pi,N));
ru=[-cos(pi*(0:N-1)/(N-1))]; %GL points
zu=[-1,-cos(pi*((0:N-2)+0.5)/(N-1)),1]; %G points + border
% y2=-cos((2*(1:N)-1)*pi/2/N)
% y2=-cos((2*(1:(N-1))-1)*pi/2/(N-1))

% clf;
% plot(xu,xu,'-x'); hold on;
% plot(yu,yu,'-o')
% plot(y2,y2,'-<')

rw=zu;
zw=ru;

rv=rw;
zv=zu;

rp=rv(2:end-1);
zp=zv(2:end-1);

%
u=zeros((N+1)^2*4,1); % as always p u v w
% J=sparse((N+1)^2*4,(N+1)^2*4,0);
% J=sparse((N+1)^2*4,(N+1)^2*4,0);
g=u; % conti umom vmom wmom

%
nsInterp;
u(3:4:end)=av;
u(2:4:end)=au;
u(1:4:end)=ap;
u(4:4:end)=aw;
%%
dlmwrite("u.dat",u,'precision',17);
system(" cp u.dat /home/gesla/Documents/git/rotst2/build/u-21.dat");


save('input');
%%

tic;
% [g,jac]=evalRhsAndJac(rp,zp,ru,zu,rv,zv,rw,zw,N,R,H,u); 
evalRhsAndJac(); 
% !time ./run_evalRhsAndJac.sh /home/gesla/mtb-installation-folder/
toc;

%
a=load("output.mat");
g=a.g;
%
% mesh(reshape(g(1:4:end),[N+1,N+1])) % conti
mesh(reshape(g(1:4:end),[N+1,N+1]))
%%
mesh(reshape(log10(abs(av)),[N+1,N+1]))


%%
% tiledlayout(1,2);
% nexttile;

cd /home/gesla/Documents/git/rotst2/build
a=importdata("g.dat"); n=12;
mesh(rv,zv,reshape(a(2:4:end),[n,n]))

%%
nexttile;
cd /home/gesla/Documents/git/mtb-work/cheb-diff/rotst
a=importdata("rhs-nu.dat"); n=153;

mesh(reshape(a(2:4:end),[n,n]))

%%
clc; clear;
% set(gcf,"Position",[   515    93   959   689])
set(gcf,"Position",[   109          93        1365         689])
varr=[];
narr=[];
%
%  cd /home/gesla/Documents/git/rotst2/build/;
 %
% N=14; 
% name="";


% cas="3"; name="-regsh"; regname=", shroud only regularisation eps=0.03"; for N=4:2:24
% cas="4"; name="-regst"; regname=", stator only regularisation eps=0.03"; for N=4:2:24

% cas="5";  name="-dblreg"; regname=", shroud + stator regularisation eps=0.03"; for N=4:2:26

cas="6";  name="-dblreg2"; regname=",  shroud + stator regularisation eps=0.006"; for N=4:2:24
% cas="1"; name="-constsh"; regname=", no regularisation"; for N=4:2:24
% cas="2"; name="-lin"; regname=""; for N=4:2:28
% for N=12
a=importdata("v-"+num2str(N)+".dat"+name);
v=reshape(a(:,3),[N+1,N+1]);

varr=[varr;v((N+2)/2,(N+2)/2)];
narr=[narr;N];
end
    clf;
tiledlayout("flow");
%(2,3);
nexttile;
a=importdata("p-"+num2str(N)+".dat"+name);
x=reshape(a(:,1),[N-1,N-1]);
y=reshape(a(:,2),[N-1,N-1]);
p=reshape(a(:,3),[N-1,N-1]);
mesh(x,y,p);
title("p"); 

%
nexttile;
a=importdata("u-"+num2str(N)+".dat"+name);
x=reshape(a(:,1),[N+1,N]);
y=reshape(a(:,2),[N+1,N]);
u=reshape(a(:,3),[N+1,N]);
mesh(x,y,u);
title("u"); 


nexttile;
a=importdata("v-"+num2str(N)+".dat"+name);
x=reshape(a(:,1),[N+1,N+1]);
y=reshape(a(:,2),[N+1,N+1]);
v=reshape(a(:,3),[N+1,N+1]);
mesh(x,y,v);
title("v"); 


nexttile;
a=importdata("w-"+num2str(N)+".dat"+name);
x=reshape(a(:,1),[N,N+1]);
y=reshape(a(:,2),[N,N+1]);
w=reshape(a(:,3),[N,N+1]);
mesh(x,y,w);
title("w"); 

% ru=[-cos(pi*(0:N-1)/(N-1))]; %GL points

nexttile;
loglog(narr,abs(varr-varr(end)),'-o'); grid on; hold on;
title("midpoint error: |v_n(0,0)-v_{max n}(0,0)|")
xlabel("n"); ylabel("error");


nexttile;
loglog(narr,abs(varr-varr(end)),'-o'); grid on; hold on;

for is=2:2:16
%     loglog(narr,narr.^(-is),'k--');
    loglog(narr,abs(varr(1)-varr(end))*(narr./narr(1)).^(-is),'k--');
%     text(narr(end),narr(end).^(-is),"slope "+num2str(-is))
    text(narr(end),abs(varr(1)-varr(end))*(narr(end)./narr(1)).^(-is),"slope "+num2str(-is))
end
title("|v_n(0,0)-v_{max n}(0,0)|");
xlabel("n"); ylabel("error");
sgtitle("Case "+num2str(cas)+" | Chebyshev solution rotor-stator, Re=10, (R=1,H=1, square cavity at the axis), N="+num2str(N)+regname)
%
exportgraphics(gcf,"chebStudy"+"-c"+cas+name+".png")
%%
% save("lin","varr","narr");
% save("dblreg","varr","narr");
save("dblreg2","varr","narr");
%%
clf;
loglog(narr,abs(varr-varr(end)),'-o'); grid on; hold on;
a=load("lin.mat"); loglog(a.narr,abs(a.varr-a.varr(end)),'-x'); 
a=load("dblreg.mat"); loglog(a.narr,abs(a.varr-a.varr(end)),'-x'); 
a=load("dblreg2.mat"); loglog(a.narr,abs(a.varr-a.varr(end)),'-x'); 
%%
close all;
N=20;
ru=[-cos(pi*(0:N-1)/(N-1))]; %GL points
plot(ru,exp(-(ru+1)/2/0.003),'-x')
%%
a=importdata("umapped.dat");
a1=importdata("u1mapped.dat");
%%
clf;
u=reshape(a(4:4:end),[15,15]);
u1=reshape(a1(4:4:end),[13,13]);
mesh(u(1:13,1:13)-u1)
% hold on;
%%
figure();
mesh()


