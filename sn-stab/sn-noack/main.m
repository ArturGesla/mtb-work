clc; close all; clear;
%%

% Noack system
% f=[mu*u-v*g-w*u;
%     mu*v+u/g-v*w;
%     -w+u^2+v^2*g*g];

clc; clear; close all;
%%
%
%parms
% it=it+1;
% mu=0.04; r=sqrt(mu); gm=1; %gamma
% mu=0.04+it*0.04; r=sqrt(mu); gm=1; %gamma
mu=3/4; r=sqrt(mu); gm=1; %gamma
neq=3; np=100; np=np+2;

%init
t=0:2*pi/(np-1):2*pi; u=r*cos(t); v=r*sin(t)/gm; w=r^2*ones(1,length(t));
T=2*pi; x=u; y=v; z=w;

X=[x(1:end-1)',y(1:end-1)',z(1:end-1)'];
X=X+rand(size(X))*1e-2;

%

z=fft(X); nt=3;
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
[g,jac]=calculateRhsAndJac(3,nt,u,mu,gm);
u=u-jac\g';
fprintf("it: %d \t norm(rhs): %4.2e\n",i,norm(g))
end

% stability

j2=jac(1:end-1,1:end-1);
[evc,evs]=eig(full(j2)); evs=diag(evs);
%
close all;
plot(evs,'x'); grid on; hold on;
text(real(evs),imag(evs),num2str([1:length(evs)]')); 
% w maire ok z fd

evsAnal=[0;(-1-sqrt(1-8*mu))/2;(-1+sqrt(1-8*mu))/2;];
plot(real(evsAnal),imag(evsAnal),'o');
xlim([-1 0.4]);
ylim([-4 4])
exportgraphics(gcf,"spec"+num2str(it)+".png");
%%
close all;
clf;

up=u; ntp=nt;
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xpb=xp;
plot(xp)
%%
close all;
clf;

up=[evc(:,11);0]; ntp=nt;
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp); xpp=xp;
% xp=xp+xpb;  plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; 

plot(xp)
%%
up=[evc(:,12);0]; ntp=nt;
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
xp=xp+xpb;  plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; 

up=[evc(:,13);0]; ntp=nt;
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
xp=xp+xpb; plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; 

% up=[evc(:,3);0]; ntp=nt2;
% zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
% xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; xpp=xp;

%% not trzeba to przemyslec
