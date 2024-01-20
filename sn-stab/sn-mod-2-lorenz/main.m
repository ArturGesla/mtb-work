clc; close all; clear;
% cd     '/people/gesla/Documents/git/mtb-work/sn-stab/sn-mod-noack';
% 


% Lorenz system
np=52; r=24;
main_lorenz_ti
%



z=fft(X); nt=17; if(mod(nt,2)==0) error("nt even"); end
arr=[1:nt]; %arr(2:2:end)=[];
a1=arr; arr=[arr,length(z)-fliplr(arr(2:end))+2];
zcut=z*0; zcut(arr,:)=z(arr,:);
X2=ifft(zcut);

close all;
% plot(X); hold on;
plot(X2); hold on

z1=zcut(a1,:)./(np-1)/2;
om=2*pi/T;
% phase
halfsol = @(phi) angle(sum([0:nt-1]'.*z1(:,1).*exp(i*phi*[0:nt-1])')); %fully wrong wtf
phi=fzero(halfsol,0);
z1=z1.*exp(i*phi*[0:nt-1])';
z11=z*0; z11(1:nt,:)=z1; z11(end-nt+2:end,:)=flipud(conj(z1(2:end,:)));
X3=ifft(z11*length(z11));
set(gca,"ColorOrderIndex",1);plot(X3,'--'); hold on

u=[reshape(real(z1.'),[3*nt,1]);reshape(imag(z1.'),[3*nt,1])];
u(nt*3*2+1)=om/2;
uinit=u;
%% load
% load("../lorenz-sn/uSN-lorenz-nt4.mat"); om=u(end); u=reshape(u(1:end-1),[3,4*2]); 
% uinit=reshape(uinit(1:end-1),[3,nt*2]); uinit=uinit*0; uinit(:,1:2:nt)=u(:,1:4); uinit(:,nt+1:2:end)=u(:,5:end);

load("../lorenz-sn/uSN-lorenz-nt9.mat"); om=u(end); u=reshape(u(1:end-1),[3,9*2]); 
uinit=reshape(uinit(1:end-1),[3,nt*2]); uinit=uinit*0; uinit(:,1:2:nt)=u(:,1:9); uinit(:,nt+1:2:end)=u(:,10:end);

u=reshape(uinit,[3*nt*2,1]); u(end+1)=om/2; 
%%
for i=1:20
[g,jac]=calculateRhsAndJac(3,nt,u,r);
u=u-jac\g';
fprintf("it: %d \t norm(rhs): %4.2e\n",i,norm(g))
end

% clf; spy(jac); grid on; grid minor;
%% second stability

nt2=nt;
% l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(1,1)=1; l31=reshape(l3,[3*nt2*2,1])';
% l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(2,1)=1; l32=reshape(l3,[3*nt2*2,1])';
% l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(3,1)=1; l33=reshape(l3,[3*nt2*2,1])';

l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2; l3(1,1)=1; r31=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2; l3(2,1)=1; r32=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2; l3(3,1)=1; r33=reshape(l3,[3*nt2*2,1])';


jmod=full(jac(1:end-1,1:end-1));%-b*j2hm+b;

bmod=jmod*0;
bmod(4,:)=r31;
bmod(5,:)=r32;
bmod(6,:)=r33;

[evc,evs]=eig(full(jmod),full(bmod)); evs=diag(evs); b=abs(evs)<Inf; evs=evs(b); evc=evc(:,b);
lam=evs+1;
log(lam)/(pi/u(end))

