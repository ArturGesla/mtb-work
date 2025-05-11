clc; close all; clear;
% cd     '/people/gesla/Documents/git/mtb-work/sn-stab/sn-mod-noack';
%%

% Lorenz system
np=52; r=24;
main_lorenz_ti
%



z=fft(X)./length(X); nt=7; if(mod(nt,2)==0) error("nt even"); end
arr=[1:nt]; %arr(2:2:end)=[];
a1=arr; arr=[arr,length(z)-fliplr(arr(2:end))+2];
zcut=z*0; zcut(arr,:)=z(arr,:);
X2=ifft(zcut)*length(zcut);

close all;
 plot(X); hold on; set(gca,"ColorOrderIndex",1);
plot(X2,':'); hold on
%
z1=zcut(1:nt,:);
om=2*pi/T;
% phase
halfsol = @(phi) real(1i*sum([0:nt-1]'.*z1(:,1).*exp(1i*phi*[0:nt-1]).')); %fully wrong wtf
phi=fzero(halfsol,0);
z1=z1.*exp(i*phi*[0:nt-1]).';
z11=z*0; z11(1:nt,:)=z1; z11(end-nt+2:end,:)=flipud(conj(z1(2:end,:)));
X3=ifft(z11)*length(z11);
set(gca,"ColorOrderIndex",1);plot(X3,'--'); 
%
u=[reshape(real(z1.'),[3*nt,1]);reshape(imag(z1.'),[3*nt,1])];
u(nt*3*2+1)=om/2;
uinit=u;
%% load
% load("../lorenz-sn/uSN-lorenz-nt4.mat"); om=u(end); u=reshape(u(1:end-1),[3,4*2]); 
% 
% uinit=reshape(uinit(1:end-1),[3,nt*2]); uinit=uinit*0; uinit(:,1:2:nt)=u(:,1:4); uinit(:,nt+1:2:end)=u(:,5:end);
% u=reshape(uinit,[3*nt*2,1]); u(end+1)=om/2; 
%%
load unt7.mat;
for i=1:7
[g,jac]=calculateRhsAndJac(3,nt,u,r);
u=u-jac\g';
fprintf("it: %d \t norm(rhs): %4.2e\n",i,norm(g))
end

% clf; spy(jac); grid on; grid minor;
% svds(jac(1:end-1,1:end-1),5,'smallest')'
%% second stability

nt2=nt;
l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(1,1)=1; l31=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(2,1)=1; l32=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(3,1)=1; l33=reshape(l3,[3*nt2*2,1])';

l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2; l3(1,1)=1; r31=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2; l3(2,1)=1; r32=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2; l3(3,1)=1; r33=reshape(l3,[3*nt2*2,1])';


jmod=full(jac(1:end-1,1:end-1));%-b*j2hm+b;
ind=nt*3+1;
jmod(ind,:)=l31;
jmod(ind+1,:)=l32;
jmod(ind+2,:)=l33;

bmod=jmod*0;
bmod(ind,:)=r31;
bmod(ind+1,:)=r32;
bmod(ind+2,:)=r33;

[evc,evs]=eig(full(jmod),full(bmod)); evs=diag(evs); b=abs(evs)<Inf; evs=evs(b); evc=evc(:,b);

fprintf("Result || \n");
fprintf("Numerical fl mult:\n");
disp(evs');

%%

up=u; ntp=nt; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
nz=(np-1)*2-length(zp)*2+1; zp=[zp;zeros(nz,3);conj(flipud(zp(2:end,:)))]; xp=ifft(zp)*(length(zp)); xp=real(xp);

clf; plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xpb=xp; plot3(xp(1,1),xp(1,2),xp(1,3),'o');plot3(xp(2,1),xp(2,2),xp(2,3),'>');
% plot(xp)

%%
xp=X; plot3(xp(:,1),xp(:,2),xp(:,3));