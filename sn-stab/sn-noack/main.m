clc; close all; clear;
%%

% Noack system
% f=[mu*u-v*g-w*u;
%     mu*v+u/g-v*w;
%     -w+u^2+v^2*g*g];

clc; clear; close all;
%
%parms
mu=0.04; r=sqrt(mu); gm=1; %gamma
neq=3; np=100; np=np+2;

%init
t=0:2*pi/(np-1):2*pi; u=r*cos(t); v=r*sin(t)/gm; w=r^2*ones(1,length(t));
T=2*pi; x=u; y=v; z=w;

X=[x(1:end-1)',y(1:end-1)',z(1:end-1)'];
X=X+rand(size(X))*1e-2;

%%

z=fft(X); nt=31;
arr=[1:nt]; a1=arr; arr=[arr,length(z)-fliplr(arr(1:end-1))+1];
zcut=z*0; zcut(arr,:)=z(arr,:);
X2=ifft(zcut);

close all;
plot(X); hold on;
plot(X2); hold on;
%
z1=zcut(a1,:)./(np-1);
om=2*pi/T;
an=angle(sum(z1(2,1))); %this could be better
z1=z1.*exp(-an*1i*[0:nt-1]');
angle(sum(z1(:,1)))
% zcut(arr,:)=zcut(arr,:).*exp(an*1i*[arr-1]');
% X3=ifft(zcut);
% plot(X3); hold on;
% x=real(z1(1,1)+z1(2,1)*exp(i*om*t)+z1(3,1)*exp(2*i*om*t)+z1(2,1)'*exp(-i*om*t)+z1(3,1)'*exp(-2*i*om*t));
x=real(z1(1,1)+z1(2,1)*exp(1i*om*t)+z1(2,1)'*exp(-1i*om*t));
plot(x,'-x')
%
% z1=zcut(a1,:);
u=[reshape(real(z1.'),[3*nt,1]);reshape(imag(z1.'),[3*nt,1])];
u(nt*3*2+1)=om;
u0=u;
%
%
for i=1:15
[g,jac]=calculateRhsAndJac(3,nt,u);
u=u-jac\g';
fprintf("it: %d \t norm(rhs): %4.2e\n",i,norm(g))
end

%% stability

j2=jac(1:end-1,1:end-1);
[evc,evs]=eig(full(j2)); evs=diag(evs);
%%
close all;
plot(evs,'x'); grid on; 
text(real(evs),imag(evs),num2str([1:length(evs)]')); 
% w maire ok z fd

%% new method - remap to half mode
nt2=(nt-1)*2+1;
u2=[]; u2(nt2*3*2+1)=om/2;

z=reshape(u(1:(end-1)/2),[3,nt])'+1i*reshape(u((end-1)/2+1:(end-1)),[3,nt])';
z2=reshape(u2(1:(end-1)/2),[3,nt2])'+1i*reshape(u2((end-1)/2+1:(end-1)),[3,nt2])';

z2(1:2:end,:)=z;

u2(1:end-1)=[real(reshape(z2.',[nt2*3,1])); imag(reshape(z2.',[nt2*3,1]))]; u2=u2.';
% u2 should solve g

[g,jac]=calculateRhsAndJac(3,nt2,u2); norm(g)

j2hm=jac(1:end-1,1:end-1);
[evc2,evs2]=eig(full(j2hm)); evs2=diag(evs2);

%%
clf;
plot(evs,'x'); grid on;
hold on; plot(evs2,'+'); grid on;
% text(real(evs2),imag(evs2),num2str([1:length(evs2)]')); 

%%
ll=mod(floor((0:nt2*3-1)./3),2);
ll2=[ll,ll];
% b=zeros(length(evs2)); b(4:6,4:6)=eye(3);
b=diag(ll2);
%
% 
% [evc3,evs3]=eig(full(j2hm),b); evs3=diag(evs3); plot(evs3,'sq'); grid on;
% [evc3,evs3]=eig(full(j2hm-b*j2hm+b),b); evs3=diag(evs3); plot(evs3,'sq'); grid on;



%
% spy(j2hm-b*j2hm+b);
% spy(b)
% 
% grid on; set(gca,"XTick",[1:3:nt2*3*2]);set(gca,"YTick",[1:3:nt2*3*2])

%
l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(1,1)=1; l31=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(2,1)=1; l32=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(3,1)=1; l33=reshape(l3,[3*nt2*2,1])';

l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2; l3(1,1)=1; r31=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2; l3(2,1)=1; r32=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2; l3(3,1)=1; r33=reshape(l3,[3*nt2*2,1])';

%
jmod=j2hm;%-b*j2hm+b;
jmod(4,:)=l31;
jmod(5,:)=l32;
jmod(6,:)=l33;

bmod=jmod*0;
bmod(4,:)=r31;
bmod(5,:)=r32;
bmod(6,:)=r33;

%
[evc,evs]=eig(full(jmod),full(bmod)); evs=diag(evs); b=abs(evs)<Inf; evs=evs(b); evc=evc(:,b);

%%
% close all;
clf;

up=u2; ntp=nt2;
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xpb=xp;

up=[evc(:,3);0]; ntp=nt2;
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; xpp=xp;

%% not trzeba to przemyslec
