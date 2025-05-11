clc; clear; uM=[]; nM=[];
%%
n=101;
nu=1/10;
dt=0.01/2/2;
tmax=1;
% nt=10;
nt=tmax/dt;
dx=1/(n-1);
x=0:dx:1;

um1=ones(n,1);
u=um1;
u=(1-cos(x'*2*pi))/2;
u=(sin(x'*2*pi));

D2=(diag(ones(n-1,1),-1)-2*eye(n)+diag(ones(n-1,1),1))/dx/dx;
A=eye(n)/dt-nu*D2;
A(1,1:2)=[1,0]; A(end,end-1:end)=[0,1];

% hold on;
for i=1:nt
% rhs=(4*u-um1)/2/dt; rhs(1)=0; rhs(end)=0;
rhs=(eye(n)/dt)*u; rhs(1)=0; rhs(end)=0;
um1=u;
u=A\rhs;
% plot(u)

end
uM=[uM,u];
nM=[nM, norm(u)];
%%
uex=exp(-tmax*nu*4*pi*pi)*sin(2*pi*x);
hold on;
plot(x,u)
plot(x,uex)
%%
ii=6;
(nM(2)-nM(1))/(nM(3)-nM(2))
(uM(ii,2)-uM(ii,1))/(uM(ii,3)-uM(ii,2))
