clc; close all; clear;
%%
np=100; r=24;
main_lorenz_ti
%%
z=fft(X); nt=4;
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

for i=1:8
[g,jac]=calculateRhsAndJac(3,nt,u);
u=u-jac\g';
norm(g)
end

%%
z2=reshape(u(1:nt*3),[3,nt])'+reshape(u(nt*3+1:end-1),[3,nt])'*1i;
% X4=ifft([z2;z2*0;flipud(z2(2:end,:))])

%%
evsM={}; leg=[];
%%
rank(full(jac(1:end-1,1:end-1)))
j2=(full(jac(1:end-1,1:end-1)));
close all;
[evc,evs]=eig(j2); evs=diag(evs);
evsM{end+1}=evs; leg=[leg;"nt = "+num2str(nt)];
plot(evs,'x'); grid on;

%%
close all;
for i=1:length(evsM)
    plot(evsM{i},'x'); grid on; hold on;
end
xlim([-0.1 0.1]);
legend(leg);
xlabel("real"); ylabel("imag");
title("spectrum of cut jacobian spectral newton Lorenz ")
%%
exportgraphics(gcf,"spec2.png");