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
% mu=0.04; r=sqrt(mu); gm=1; %gamma
% mu=0.04+it*0.04; r=sqrt(mu); gm=1; %gamma
mu=3/4; r=sqrt(mu); gm=1; %gamma
neq=3; np=100; np=np+1;

%init
t=0:4*pi/(np-1):2*pi*2; u=r*cos(t); v=r*sin(t)/gm; w=r^2*ones(1,length(t));
T=2*pi; x=u; y=v; z=w;

X=[x(1:end-1)',y(1:end-1)',z(1:end-1)'];
% X=X+rand(size(X))*1e-2;

%

z=fft(X); nt=11; if(mod(nt,2)==0) error("nt even"); end
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

% %% stability
% 
% % u((end-1)/2+1:end)=-u((end-1)/2+1:end);
% % u((end-1)/2+1:end-1)=-u((end-1)/2+1:end-1);
% % [g,jac]=calculateRhsAndJac(3,nt,u,mu,gm);
% 
% j2=jac(1:end-1,1:end-1);
% [evc,evs]=eig(full(j2)); evs=diag(evs);
% %
% close all;
% plot(evs,'x'); grid on; hold on;
% text(real(evs),imag(evs),num2str([1:length(evs)]')); 
% % w maire ok z fd
% 
% evsAnal=[0;(-1-sqrt(1-8*mu))/2;(-1+sqrt(1-8*mu))/2;];
% flmult=exp(evsAnal*T);
% plot(real(evsAnal),imag(evsAnal),'o');
% xlim([-1 0.4]);
% ylim([-4 4])
% % exportgraphics(gcf,"spec1.png");
% %%
% close all;
% clf;
% %%
% up=u; ntp=nt; 
% zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
% xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xpb=xp;
% 
% plot3(xp(1,1),xp(1,2),xp(1,3),'o');plot3(xp(2,1),xp(2,2),xp(2,3),'>');
% % plot(xp)
% %%
% close all;
% clf;
% iev=20;
% evs(iev)
% 
% % up=[evc(:,iev);0]; ntp=nt; u1=up(1:end-1);
% % up=[j2*u2;0]; ntp=nt; 
% up=[u4;0]; ntp=nt; 
% zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
% xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp); xpp=xp;
% 
% xp=xpp+xpb;  
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; 
% xp=xpb;  
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; 
% 
% % plot(xpp+xpb)
% 
% %% expo counterpart
% ll=length(xp); sigma=evs(iev);
% t=0:2*pi/ll:2*pi-2*pi/ll;
% xp=xpp.*(exp(sigma*t)');
% 
% z2=fft(xp)./length(xp);
% u2=z2(1:nt,:); 
% u2=[reshape(real(u2.'),[3*nt,1]);reshape(imag(u2.'),[3*nt,1])]; 
% %%
% clf; ip=1; plot(up(ip:3:end-1),'x-'); hold on; plot(u2(ip:3:end),'x-');  u3=j2*u2; plot(u3(ip:3:end),'x-')
% %%
% up=[evc(:,12);0]; ntp=nt;
% zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
% xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
% xp=xp+xpb;  plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; 
% 
% up=[evc(:,13);0]; ntp=nt;
% zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
% xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
% xp=xp+xpb; plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; 
% 
% % up=[evc(:,3);0]; ntp=nt2;
% % zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
% % xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
% % plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; xpp=xp;

%% not trzeba to przemyslec
