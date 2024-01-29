clc; close all; clear;
% cd     '/people/gesla/Documents/git/mtb-work/sn-stab/sn-mod-noack';
%%

% Noack system
% f=[mu*u-v*g-w*u;
%     mu*v+u/g-v*w;
%     -w+u^2+v^2*g*g];

% clc; clear; close all;
%
%
%parms
% it=it+1;
mu=0.04; r=sqrt(mu); gm=1; %gamma
% mu=0.04+it*0.04; r=sqrt(mu); gm=1; %gamma
% mu=3/4; r=sqrt(mu); gm=1; %gamma
neq=3; np=100; np=np+1;

%init
t=0:4*pi/(np-1):2*pi*2; u=r*cos(t); v=r*sin(t)/gm; w=r^2*ones(1,length(t));
T=2*pi; x=u; y=v; z=w;

X=[x(1:end-1)',y(1:end-1)',z(1:end-1)'];
% X=X+rand(size(X))*1e-2;

%

z=fft(X); nt=5; if(mod(nt,2)==0) error("nt even"); end
arr=[1:nt]; a1=arr; arr=[arr,length(z)-fliplr(arr(1:end-1))+1];
zcut=z*0; zcut(arr,:)=z(arr,:);
X2=ifft(zcut);

close all;
plot(X); hold on;
plot(X2); hold on;
%
z1=zcut(a1,:)./(np-1);
om=2*pi/T;
% an=angle(sum(z1(2,1))); %this could be better
% z1=z1.*exp(-an*1i*[0:nt-1]');
% angle(sum(z1(:,1)))
% zcut(arr,:)=zcut(arr,:).*exp(an*1i*[arr-1]');
% X3=ifft(zcut);
% plot(X3); hold on;
% x=real(z1(1,1)+z1(2,1)*exp(i*om*t)+z1(3,1)*exp(2*i*om*t)+z1(2,1)'*exp(-i*om*t)+z1(3,1)'*exp(-2*i*om*t));
% x=real(z1(1,1)+z1(2,1)*exp(1i*om*t)+z1(2,1)'*exp(-1i*om*t));
plot(x,'-x')
%
% z1=zcut(a1,:);
u=[reshape(real(z1.'),[3*nt,1]);reshape(imag(z1.'),[3*nt,1])];
u(nt*3*2+1)=om/2;
% u0=u; u(end)=1
% load("../sn-noack/u2.mat"); u=u2;
% load("../sn-noack/u2-nt3.mat"); u=u2;
%
%
for i=1:15
[g,jac]=calculateRhsAndJac(3,nt,u,mu,gm);
u=u-jac\g';
fprintf("it: %d \t norm(rhs): %4.2e\n",i,norm(g))
end

% clf; spy(jac); grid on; grid minor;
svds(jac(1:end-1,1:end-1),5,'smallest')
% second stability

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

evsAnal=[0;(-1-sqrt(1-8*mu))/2;(-1+sqrt(1-8*mu))/2;];
flmult=exp(evsAnal*T);

fprintf("Result || \n");
fprintf("Analytical fl mult:\n");
disp(flmult');
fprintf("Numerical fl mult:\n");
disp(evs');

%% visu
% up=u; ntp=nt; 
up=u+[evc(:,3);0]; ntp=nt; 
% up=[evc(:,1);0]; ntp=nt; 
% up=u+[V(:,5)]; ntp=nt; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
nz=(np-1)*2-length(zp)*2+1; zp=[zp;zeros(nz,3);conj(flipud(zp(2:end,:)))]; xp=ifft(zp)*(length(zp)); xp=real(xp);

clf; plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xpb=xp; plot3(xp(1,1),xp(1,2),xp(1,3),'o');plot3(xp(2,1),xp(2,2),xp(2,3),'>');
% clf; plot(xp)