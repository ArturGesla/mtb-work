clf; leg=[];

a=load('eva-omline-k-1-l-120.mat'); plot(a.eva,'k.'); hold on; leg=[leg;max(a.x)];
a=load('eva-omline-k-1-l-130.mat'); plot(a.eva,'b.'); hold on; leg=[leg;max(a.x)];
a=load('eva-omline-k-1-l-140.mat'); plot(a.eva,'r.'); hold on; leg=[leg;max(a.x)];

legend("zmax="+num2str(leg));

xlim([0 1]); ylim([-0.6 0.6 ]); xlabel("$a_r$"); ylabel("$a_i$"); grid on;
title("Re="+num2str(R)+", $\beta$="+num2str(bbar)+", $\omega\in$("+num2str(min(omegaa))+", "+num2str(max(omegaa))+")");
size_sq23; fnts=10; jfm_plt_aid_comm;

exportgraphics(gcf,"p4-lw97fig6-comp1.eps")%lw97 fig 6

%%

clf;leg=[];

a=load('eva-areal-k-1-l-120+.mat'); plot(a.eva,'k.'); hold on; leg=[leg;max(a.x)];
a=load('eva-areal-k-1-l-130+.mat'); plot(a.eva,'b.'); hold on; leg=[leg;max(a.x)];
a=load('eva-areal-k-1-l-140+.mat'); plot(a.eva,'r.'); hold on; leg=[leg;max(a.x)];

legend("zmax="+num2str(leg),"Location","best");


xlim([-0.08 0.04]); ylim([-0.08 0.01 ]); xlabel("$\omega_r$"); ylabel("$\omega_i$");  grid on;
title("Re="+num2str(R)+", $\beta$="+num2str(bbar)+", $\alpha\in$("+num2str(min(ar))+", "+num2str(max(ar))+")");
size_sq23; fnts=10; jfm_plt_aid_comm;

exportgraphics(gcf,"p4-lw97fig6-comp2.eps")%lw97 fig 6

