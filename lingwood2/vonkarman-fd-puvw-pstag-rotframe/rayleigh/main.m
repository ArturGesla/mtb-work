 clc; clear;
%%
% a=load("../vk-np-140.mat");
a=load("../vk-np-100.mat");
% a=load("../vk-np-200.mat");
% a=load("../vk-np-400.mat");
% a=load("../vk-np-800.mat");
x=a.x;
u=a.u*0;
U=a.u;
omega=0+0.04i;
% beta=bbar=beta/R;
bbar=0.126;
R=1;
alpha=1; 

[g,jac0,jac1,jac2,jac3]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);

%%
eva=[];
omega=-0.0262+0.015i;
oma=-0.0692 :0.002/4:0.0262; oma1=oma+0.04i;
oma=-0.0692 :0.002/4:0.0262; oma2=oma+0.015i;
% om0=-0.0692; om1=0.0262; omx=-1:2/50:1; omx=omx.^3; omx=(omx+1)/2*(om1-om0)+om0; oma3=omx+0.0133i;
oma=-0.0692 :0.002/4:0.0262; oma3=oma+0.0133i;
%
% ev=eig(full(jac0),-full(jac1));

%%
% omega=-0.0262+0.0125i;
% omega=-0.0262+0.04i;
% omega=-0.0262+0.015i;
omega=omega-0.002;
[g,jac0,jac1,jac2,jac3]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
[evc,evs]=polyeig(jac0,jac1,jac2,jac3); ev=(evs);
% [evc,evs]=eig(full(jac0),-full(jac1)); ev=diag(evs);
% [evc,evs]=eigs((jac0),-(jac1),20,'smallestabs'); ev=diag(evs);
eva=[eva,ev];
%%
clf;
plot(eva,'k.'); hold on;
%  plot(ev,'x'); hold on;; text(real(ev),imag(ev),num2str([1:length(ev)]'))
% plot(ev,'o');
xlim([0 0.45]); ylim([-0.5 0.5]); grid on;
% xlim([-0.1 1.5]); ylim([-1 1]); grid on;
% xlim([-0.1 0.5]*10); ylim([-0.5 0.5]*10); grid on;
% xlim([-0.1 0.5]*40); ylim([-0.5 0.5]*40)

pbaspect([5 1 1]); xlabel("a_r"); ylabel("a_i"); 
% fnts=12; jfm_plt_aid_comm;
% exportgraphics(gcf,"spatialbranches-c.eps")
% exportgraphics(gcf,"spatialbranches-b.eps")
exportgraphics(gcf,"spatialbranches-a.eps")

%%
clf;
iev=375;
up=reshape(real(evc(:,iev)),[1,length(x)])';
plot(up(:,1:end),x,'x-');
title("ev "+num2str(iev)+":"+num2str(ev(iev),'%4.2e'))
iev=iev+1
% ylim([0 0.8])
%%
eva=[];
oma=oma1;
for i=1:length(oma)
%     for i=40:length(oma)
omega=oma(i);
[g,jac0,jac1,jac2,jac3]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
[evc,evs]=polyeig(jac0,jac1,jac2,jac3); ev=(evs);
% [evc,evs]=eig(full(jac0),-full(jac1)); ev=diag(evs);
% [evc,evs]=eigs((jac0),-(jac1),20,'smallestabs'); ev=diag(evs);
eva=[eva,ev];
i
end