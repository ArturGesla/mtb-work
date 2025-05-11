clc; clear;
% %%
% lam=-0.1:0.001:0;
% b0=0.3;
% gamma=-10;
% beta=pi/4;
% 
% 
% t1=real((acos((b0*cos(beta)-lam)/(b0))+beta+2*pi)./(1-b0*sin(beta)+sqrt(lam.*(2*b0*cos(beta)-lam)+b0^2*sin(beta)^2))/pi);
% t2=real((acos((b0*cos(beta)-lam)/(b0))+beta+2*pi)./(1-b0*sin(beta)-sqrt(lam.*(2*b0*cos(beta)-lam)+b0^2*sin(beta)^2))/pi);
% t3=real((-acos((b0*cos(beta)-lam)/(b0))+beta+2*pi)./(1-b0*sin(beta)+sqrt(lam.*(2*b0*cos(beta)-lam)+b0^2*sin(beta)^2))/pi);
% t4=real((-acos((b0*cos(beta)-lam)/(b0))+beta+2*pi)./(1-b0*sin(beta)-sqrt(lam.*(2*b0*cos(beta)-lam)+b0^2*sin(beta)^2))/pi);
% 
% phi=0:0.01:20;
% l1=b0*(cos(beta)-cos(beta-phi));
% l2=phi./(1+b0*(sin(beta-phi)-sin(beta)));
% 
% l3=2./(1-gamma*lam);
% %%
% hold on;
% plot(lam,t1)
% plot(lam,t2)
% plot(lam,t3)
% plot(lam,t4)
% grid on; 
% plot(l1,l2/pi)
% plot(lam,l3)
% ylim([0 6])
% % xlim([-0.1 0.1])
%
% figure()
hold on;
grid on; grid minor;

phi=0:0.01:10; beta=pi/4; gamma=-10;%kc;
kc=-1/(2*pi*(gamma*sin(beta)+cos(beta)));
b0=0.02;


l1=b0*(cos(beta)-cos(beta-phi)); %lam
l2=phi./(1+b0*(sin(beta-phi)-sin(beta))); %tau_hopf
l3=pi*2./(1-gamma*l1); %tau_pyr
l4=pi*2./(1+gamma*sqrt((-1+sqrt(1+4*l1))/-2).^2);



crit=1+b0*l2.*(cos(beta-phi)+gamma*sin(beta-phi));
crit2=1+b0*l3.*(cos(beta)+gamma*sin(beta));



plot(l1,l2/pi.*(1./(crit>0)),"LineWidth",2)
plot(l1,l2/pi.*(1./(crit<0)),'--',"LineWidth",2)
plot(l1,l3/pi.*(1./(crit2>0)),".")
plot(l1,l3/pi.*(1./(crit2<0)),'.')

legend("sub","super","tc unst","tc st")

% phi=0:0.01:10; beta=pi/2; gamma=-1;%kc;
% kc=-1/(2*pi*(gamma*sin(beta)+cos(beta)));
% b0=0.3;
% 
% l1=b0*(cos(beta)-cos(beta-phi));
% l2=phi./(1+b0*(sin(beta-phi)-sin(beta)));
% crit=1+b0*l2.*(cos(beta-phi)+gamma*sin(beta-phi));
% 
% l3=2./(1-gamma*l1);
% l4=2./(1+gamma*sqrt((-1+sqrt(1+4*l1))/-2).^2);
% 
% plot(l1,l2/pi.*(1./(crit>0)))
% plot(l1,l2/pi.*(1./(crit<0)),'--')

% plot(l1,l4)



title("K "+num2str(b0)+" beta = "+num2str(beta)+" gamma = "+num2str(gamma))

% phi=0:0.01:20; b0=0.016; beta=pi/4; %kc;
% crit=1+b0*l2.*(cos(beta-phi)+gamma*sin(beta-phi));
% l1=b0*(cos(beta)-cos(beta-phi));
% l2=phi./(1+b0*(sin(beta-phi)-sin(beta)));
% plot(l1,l2/pi.*(1./(crit>0)))
% plot(l1,l2/pi.*(1./(crit<0)),'--')

% phi=0:0.01:20; b0=0.4; beta=pi/4/0.5; %kc;
% crit=1+b0*l2.*(cos(beta-phi)+gamma*sin(beta-phi));
% l1=b0*(cos(beta)-cos(beta-phi));
% l2=phi./(1+b0*(sin(beta-phi)-sin(beta)));
% plot(l1,l2/pi.*(1./(crit>0)))
% plot(l1,l2/pi.*(1./(crit<0)),'--')

% phi=0:0.01:20; b0=0.15;
% l1=b0*(cos(beta)-cos(beta-phi));
% l2=phi./(1+b0*(sin(beta-phi)-sin(beta)));
% plot(l1,l2/pi)


% xlim([-0.1 0.1])