clc; clear;

% np=400; x=linspace(0,20,np); x=(x/20).^2*20;
L=30*16; np=1600; x=linspace(0,L,np); x=(x/L).^1*L;
% x=0:0.1:20; np=length(x);
% al=0.01; x=(exp(al*x)./exp(al*20)*2-1)*20;
u=zeros(np*4,1);
 u(1:4:end)=-1;
u(2:4:end)=0;
u(3:4:end)=1;
u(4:4:end)=1;
% u(4:4:end)=1;


% a=load("../vonkarman-bvp4c/vk.mat");
% u(2:4:end)=interp1(a.sol.x,a.sol.y(1,:),x);
% u(3:4:end)=interp1(a.sol.x,a.sol.y(2,:),x)-1;
% u(4:4:end)=interp1(a.sol.x,a.sol.y(3,:),x);
% 
% [g,jac]=evalJacRhs(u,x);
% 
% norm(g)
%
% u=rand(np*4,1);
k=1;
% k=0.313;

for i=1:15
    [g,jac,~,~]=evalJacRhs(u,x,k);
    fprintf("%d i \t norm(rhs): %4.2e \t rms norm: %4.2e \n",i,norm(g),norm(g));
    u=u-jac\g;
    if(norm(g)<1e-13) break; end;

end
%
save("vk-np-"+num2str(np)+"-k-"+num2str(k)+"-L-"+num2str(L)+".mat",'u','np','x','k');
    [g,jac,zc,zw]=evalJacRhs(u,x,k);

%%
% plot(x,reshape(u,[3,np])','x-')
up=reshape(u,[4,np])';

plot(up(:,1:end),x,'x-'); 
% xlim([-1.5 0.5]); 
% ylim([0 20]); pbaspect([10 3 1]);
legend("P","F","G","H");
% legend("F","G","H");
% fnts=12; jfm_plt_aid_comm;
% exportgraphics(gcf,"p7-vel.eps")
% save("vk-np-"+num2str(np),'u','np');

%%
clf;
plot(up(:,1:end-1),zc,'+-'); hold on; 
plot(up(:,end),zw,'+-'); 
legend("P","F","G","H","Location","best");
ylim([0 30]); grid on;
title("B\"+'"'+"odewadt layer solution");
ylabel("z"); xlabel("$\phi(z)$");
fnts=12; jfm_plt_aid_comm; size_sq23;
exportgraphics(gcf,"p4-bode.eps")

%%
[x(1:11)',up(1:11,2:end)]


%% convergence
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

%%

a=load("../bode-specP/vk-np-100.mat")
