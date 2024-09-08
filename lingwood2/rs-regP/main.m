clc; clear;
%
% np=400; x=linspace(0,20,np); x=(x/20).^2*20;
%L=1; np=100+L; x=linspace(0,L,np); x=(x/L).^1*L; Re=1000; s="x"; reh=Re;
 L=round(sqrt(1000)); np=100+L; x=linspace(0,L,np); x=(x/L).^1*L; Re=1; s="x"; reh=L.^2;
disp("Reh: "+num2str(reh));
% x=0:0.1:20; np=length(x);
% al=0.01; x=(exp(al*x)./exp(al*20)*2-1)*20;
u=zeros(np*4,1);
 u(1:4:end)=0;
u(2:4:end)=0;
u(3:4:end)=0.313;
u(end)=0.313;
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

%
% for ii=1:1%100
for i=1:25
    [g,jac]=evalJacRhs(u,x,Re);
    fprintf("%d i \t norm(rhs): %4.2e \t rms norm: %4.2e \n",i,norm(g),norm(g));
    u=u-jac\g;
    if(norm(g)<1e-13) break; end;

end
% end
%
% save("vk-np-"+num2str(np),'u','np','x');

%%
clf;
% plot(x,reshape(u,[3,np])','x-')
up=reshape(u,[4,np])';
% hold on;
plot(up(:,1:end),x,'x-'); 
% plot(up(:,2:end)*sqrt(Re),x./L,"sq-"); 
% xlim([-1.5 0.5]); 
% ylim([0 20]); pbaspect([10 3 1]);
% legend("P","F","G","H");
legend("P","F","G","H");
% fnts=12; jfm_plt_aid_comm;
% exportgraphics(gcf,"p7-vel.eps")
k=1;
save("rs-np-"+num2str(np)+"-k-"+num2str(k)+"-L-"+num2str(L)+".mat",'u','np','x','k');


%% comp Re

clf;
% plot(x,reshape(u,[3,np])','x-')
up=reshape(u,[4,np])';
plot(up(:,3),x,'x-');

a=load("../bode-regP/vk-np-130-k-1.mat");
hold on; 
up=reshape(a.u,[4,length(a.u)/4])';
% plot(0.313*up(:,1:end),a.x/sqrt(Re*0.313),'k-'); 
plot(0.313*up(:,3),a.x/sqrt(Re*0.313),'k+-'); 
xlim([0,1]); ylim([0,1]); grid on;
title("Rotor-stator solution Re=1000");
ylabel("z"); xlabel("G(z)");
legend("rotor-stator","rescaled B\"+'"'+"odewadt","Location","best");
fnts=12; jfm_plt_aid_comm; size_sq23;
exportgraphics(gcf,"p4-rs-vs-bode.eps")

%% comp reh
a=load("../bode-regP/vk-np-130-k-1.mat");
hold on; 
up=reshape(a.u,[4,length(a.u)/4])';
plot(0.313*up(:,1:end),a.x/sqrt(0.313),'k-'); 
%%