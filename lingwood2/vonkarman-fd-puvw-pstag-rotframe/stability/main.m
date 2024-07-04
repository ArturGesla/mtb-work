 clc; clear;

a=load("../vk-np-100.mat");
% a=load("../vk-np-200.mat");
% a=load("../vk-np-300.mat");
% a=load("../vk-np-400.mat");
x=a.x;
u=a.u*0;
U=a.u;
omega=0+0.04i;
% beta=bbar=beta/R;
bbar=0.126;
R=1;
alpha=1; 

[g,jac0,jac1]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
eva=[];
%
% ev=eig(full(jac0),-full(jac1));

%
omega=-0.0262+0.0125i;
[g,jac0,jac1]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
[evc,evs]=eig(full(jac0),-full(jac1)); ev=diag(evs);
eva=[eva,ev];
%%
clf;
plot(ev,'x'); hold on;
text(real(ev),imag(ev),num2str([1:length(ev)]'))
% plot(ev,'o');
xlim([0 0.5]); ylim([-0.5 0.5])
% xlim([-0.1 0.5]*10); ylim([-0.5 0.5]*10)

%%
clf;
iev=385;
up=reshape(real(evc(:,iev)),[4,length(x)])';
plot(up(:,2:end),x,'x-');
title("ev "+num2str(iev)+":"+num2str(ev(iev),'%4.2e'))
iev=iev+1