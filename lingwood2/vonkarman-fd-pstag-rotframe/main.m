clc; clear;

np=220;
x=linspace(0,20,np);
x=(x/20).^1*20;
% al=0.01; x=(exp(al*x)./exp(al*20)*2-1)*20;
u=zeros(np*3,1);

a=load("../vonkarman-bvp4c/vk.mat");
u(1:3:end)=interp1(a.sol.x,a.sol.y(1,:),x);
u(2:3:end)=interp1(a.sol.x,a.sol.y(2,:),x)-1;
u(3:3:end)=interp1(a.sol.x,a.sol.y(3,:),x);

[g,jac]=evalJacRhs(u,x);

norm(g)
%%
u=rand(np*3,1);

%
for i=1:15
    [g,jac]=evalJacRhs(u,x);
    fprintf("%d i \t norm(rhs): %4.2e \t rms norm: %4.2e \n",i,norm(g),rms(g));
    u=u-jac\g;
    if(norm(g)<1e-14) break; end;

end
save("vkrf-np-"+num2str(np),'u','np','x');

%%
% plot(x,reshape(u,[3,np])','x-')
plot(reshape(u,[3,np])',x,'x-'); xlim([-1.5 0.5]); ylim([0 20]); pbaspect([10 3 1]);
legend("F","G","H");
fnts=12; jfm_plt_aid_comm;
exportgraphics(gcf,"p7-vel.eps")
% save("vk-np-"+num2str(np),'u','np');
%%
winf=[]; npa=[]; 
a=load("vk-np-50.mat"); winf=[winf;a.u(end)]; npa=[npa; a.np]
a=load("vk-np-100.mat"); winf=[winf;a.u(end)]; npa=[npa; a.np]
a=load("vk-np-200.mat"); winf=[winf;a.u(end)]; npa=[npa; a.np]
a=load("vk-np-400.mat"); winf=[winf;a.u(end)]; npa=[npa; a.np]
a=load("vk-np-800.mat"); winf=[winf;a.u(end)]; npa=[npa; a.np]
a=load("vk-np-1600.mat"); winf=[winf;a.u(end)]; npa=[npa; a.np]

orders=log((winf(1:end-2)-winf(2:end-1))./(winf(2:end-1)-winf(3:end)))/log(2)

%%
[U,S,V]=svd(full(jac));
max(diag(S))/min(diag(S))
cond(full(jac)) 
ss=diag(S);
%%
ii=1;
up=V(:,ii); %kernel
% up=U(:,ii);  %image
upp=reshape(up,[3,np])';
% semilogx(x,upp(:,1),'x-');
semilogx(x,upp,'x-');
title("sing val: "+num2str(ss(ii),'%4.2e'))
