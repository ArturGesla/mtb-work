addpath     'C:\Users\Izabela\Documents\git\rotst2\scripts\source_for_mtb'
cd     'C:\Users\Izabela\Documents\git\mtb-work\meca-pres-2024'

%% hopf ti

nP=1; r=24; np=200;
lorenz = @(t,x) [10*(x(2)-x(1));r*x(1)-1*x(1)*x(3)-x(2);x(1)*x(2)-(8/3)*x(3)]; % Anonymous Function
% [T,X] = ode45(lorenz, [0 :T/(np):1*T], [-6.2262, -11.0027,13.0515]);
ev=[ -0.0328   -0.6228    0.6818]*5;
T=0.6779; T= 0.6803*1.2; [t,X] = ode45(lorenz, [0 :T/(np-2):T/(np-2)*(np-1)*nP],[10.135982315094342  10.189521543725682  25.691556187487929]+ev);
T=0.6779; T= 0.6803*1.5; [t,X2] = ode45(lorenz, [0 :T/(np-2):T/(np-2)*(np-1)*nP],[10.135982315094342  10.189521543725682  25.691556187487929]+ev*0);
% T=0.6779; T= 0.6803; [t,X] = ode45(lorenz, [0 :T/(np-2):T/(np-2)*(np-1)*nP], [9.3 11.1 22.2]);
%
close all;
plot(X(:,1),X(:,2)); hold on;
plot(X2(:,1),X2(:,2));
xlabel("x"); ylabel("y"); grid on;

fnts=10; jfm_plt_aid_comm; size_sq23;
%%
% exportgraphics(gcf,"hopforb24.eps")
exportgraphics(gcf,"hopforb24-shoot.eps")

%%
a=load("udat.mat"); np=(length(a.u)-1)/3;
X=reshape(a.u(1:end-1),[3,np])';

close all;
plot(X(:,1),X(:,2)); hold on;
plot(X(:,1),X(:,2),'x'); hold on;
% plot(X2(:,1),X2(:,2));
xlabel("x"); ylabel("y"); grid on;

fnts=10; jfm_plt_aid_comm; size_sq23;
%%
exportgraphics(gcf,"hopforb24-coll.eps")