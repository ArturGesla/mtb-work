clc; clear;
uarr=[];
%
nx=5; 
narr=[];
%%
% nx=(nx-1)*2+1;
ny=nx;
x=linspace(-1,1,nx)'; y=x;
dx=mean(diff(x));
x=-dx/2:dx:1+dx/2;
%%

%d2udx2-1=0
g=zeros(nx*ny,1); %laplacian d2x+d2y
u=g; %u d2udx2+d2udy2
J=sparse(nx*ny,nx*ny,0);
%
% f=x*0-1*x.^2+1;
% f=x*0-1*x+1;
f=u*0;
% f=cos(pi*x);
%

g=zeros(nx*ny,1);
% u=g-1;
J=sparse(nx*ny,nx*ny,0);
B=J;


%
tic;
% u=u-J\g;
u=J\(B*f);
toc;
% norm(g)
%
%
uPhys=reshape(u,[ny,nx]);
mesh(x,y,uPhys);
%
ix=(length(x)+1)/2; iy=(length(y)+1)/2; ip=iy+(ix-1)*ny; 
 uarr=[uarr;u(ip)]
 narr=[narr;nx];
 nx=(nx-1)*2+1;
% anal=-1/pi/pi*2
%%
save('centerFD2','narr','uarr');

%%
% clf; 
B=zeros(nx,ny); B(2:end-1,2:end-1)=1; B=reshape(B,[nx*ny,1]); B=sparse(1:length(B),1:length(B),B);
% [ev,evs]=eig(full(J)); evs=diag(evs); 
[ev,evs]=eigs((J),B,20,2); evs=diag(evs);
save('evsFD','evs')
% plot(evs,'x');

uarr=[uarr;real(evs(1:5)')]
narr=[narr;nx];
%%
save('evsFD',"uarr","narr");

%%
% i=1; 
u=real(ev(:,i));
uPhys=reshape(u,[ny,nx]);
mesh(x,y,uPhys);
i=i+1;
%%
d=diff(uarr)
d(1)/d(2)
%%
clf;
loglog(narr,abs(uarr-2/pi/pi),'x-'); hold on;
loglog(narr,narr.^(-2),'-')
loglog(narr,narr.^(-4),'-')
