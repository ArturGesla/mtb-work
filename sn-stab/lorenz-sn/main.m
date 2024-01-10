clc; close all; clear;
%%
np=100; r=24;
main_lorenz_ti
%
z=fft(X); nt=7;
arr=[1:nt]; a1=arr; arr=[arr,length(z)-fliplr(arr(1:end-1))+1];
zcut=z*0; zcut(arr,:)=z(arr,:);
X2=ifft(zcut);
%
% close all;
% plot(X); hold on;
% plot(X2); hold on;
%
z1=zcut(a1,:)./np;
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

for i=1:15
[g,jac]=calculateRhsAndJac(3,nt,u);
u=u-jac\g';
norm(g)
end

%%
% close all;
clf;
z2=reshape(u(1:nt*3),[3,nt])'+reshape(u(nt*3+1:end-1),[3,nt])'*1i;
nnt=100; X4=ifft([z2;zeros(nnt,3);conj(flipud(z2(2:end,:)))]); X4=real(X4)*length(X4);
% plot3(X4(:,1),X4(:,2),X4(:,3)); hold on; plot3(X(:,1),X(:,2),X(:,3));
plot(X4(:,1),X4(:,2)); hold on; plot(X(:,1),X(:,2));

%%
close all;
% evsM={}; leg=[];
rank(full(jac(1:end-1,1:end-1)))
j2=(full(jac(1:end-1,1:end-1)));
    B=j2*0; 
for i=[1:nt*3,nt*3+4:nt*3*2]
    B(i,i)=1; 
end
% 
[evc,evs]=eig(j2,B); evs=diag(evs);
% evsM{end+1}=evs; leg=[leg;"nt = "+num2str(nt)];

om=2*pi/T;
om=u(end);
 plot(real(evs),imag(evs)/om,'+'); grid on; hold on;
 text(real(evs),imag(evs)/om,num2str([1:length(evs)]'));
 %% new idea
 figure(2); plot(real(evc(1:3:end,[21,30,34,38,35,31,22])))
%%
close all;
sym=['x','o','>','<','v','^'];
for i=1:length(evsM)
%     plot(real(evsM{i}),imag(evsM{i})/om,sym(i)); grid on; hold on;
    lam=exp(T*evsM{i});
    plot(real(lam),imag(lam),sym(i)); grid on; hold on;
end
xlim([-1.2 1.2]); axis square; axis equal; plot(cos([0:0.01:2*pi]),sin([0:0.01:2*pi]));
legend(leg);
xlabel("real"); ylabel("imag");
title("spectrum of cut jacobian spectral newton Lorenz ")
%%
exportgraphics(gcf,"spec2.png");
%% visu
close all;
f=figure("Position",[400,350,800,800*0.4]);

subplot(1,3,1);
om=u(end);  T=2*pi/om;  neq=3; t=linspace(0,T,100); x=0; y=0; z=0; %nt=3;
iev=32; uev=evc(:,iev); f=0; g=0; h=0; eev=((evs(iev)));% eev=1i*imag((evs(iev)));
% iev2=2; uev2=evc(:,iev2); f2=0; g2=0; h2=0; eev2=(evs(iev2));

% uev=(uev+uev2)/2;
% eev=(eev+eev2)/2;

for ik=0:nt-1
    mn=(ik==0)*1+1;
    
    coeff=(u(ik*neq+1)+1i*u(nt*neq+ik*neq+1))/mn;
    x=x+coeff*exp(ik*1i*om*t)+coeff'*exp(ik*-1i*om*t);
    
    coeff=(uev(ik*neq+1)+1i*uev(nt*neq+ik*neq+1))/mn;
    f=f+(coeff*exp(ik*1i*om*t)+coeff'*exp(ik*-1i*om*t)).*exp(eev.*t);
    
%     coeff=(uev2(ik*neq+1)+1i*uev2(nt*neq+ik*neq+1))/mn;
%     f2=f2+(coeff*exp(ik*1i*om*t)+coeff'*exp(ik*-1i*om*t)).*exp(eev2.*t);
    
%     coeff=(uev(ik*neq+1)+1i*uev(nt*neq+ik*neq+1))/mn;
%     f=f+(coeff'*exp(ik*1i*2*pi*t)+coeff*exp(ik*-1i*2*pi*t)).*exp(evs(iev)'*T.*t);
    
    coeff=(u(ik*neq+2)+1i*u(nt*neq+ik*neq+2))/mn;
    y=y+coeff*exp(ik*1i*om*t)+coeff'*exp(ik*-1i*om*t);
    
    coeff=(uev(ik*neq+2)+1i*uev(nt*neq+ik*neq+2))/mn;
    g=g+(coeff*exp(ik*1i*om*t)+coeff'*exp(ik*-1i*om*t)).*exp(eev.*t);
    
    coeff=(u(ik*neq+3)+1i*u(nt*neq+ik*neq+3))/mn;
    z=z+coeff*exp(ik*1i*om*t)+coeff'*exp(ik*-1i*om*t);
    
    coeff=(uev(ik*neq+3)+1i*uev(nt*neq+ik*neq+3))/mn;
    h=h+(coeff*exp(ik*1i*om*t)+coeff'*exp(ik*-1i*om*t)).*exp(eev.*t);
    
    
%     x=x+mn*u(ik*neq+1)*cos(ik*2*pi*t);%+u(nt*neq+ik*neq+1)*1i*sin(ik*2*pi*t);
%     y=y+mn*u(ik*neq+2)*cos(ik*2*pi*t);%+u(nt*neq+ik*neq+2)*1i*sin(ik*2*pi*t);
end
plot(X(:,1),X(:,2)); hold on;
plot(x+real(f),y)
plot(x,y); grid on;


subplot(1,3,2)

plot3(x,y,z);
grid on; hold on;
amp=1;
% plot3(x+amp*f,y+amp*g,z+amp*h);
plot3(x+amp*real(f),y+amp*real(g),z+amp*real(h));
% plot3(x+amp*f,y,z);