

clf; set(gcf,"Position",[          450          295         1035       346.33]);

subplot(1,4,1:2)
hold on;
a=load("eva-areal-k-1-np-100-L-20.mat"); plot(reshape(a.eva,[30*201,1]),'.');
a=load("eva-areal-k-1-np-150-L-30.mat"); plot(reshape(a.eva,[30*201,1]),'.');
a=load("eva-areal-k-1-np-200-L-40.mat"); plot(reshape(a.eva,[30*201,1]),'.');


lg="$z_{max}$ = "+num2str([20:10:40]');



xlim([-0.08 0.04]); ylim([-0.08 0.01]); grid on;

xlabel("$\omega_r$"); ylabel("$\omega_i$"); 
title("Re="+num2str(27.4)+", $\beta$="+num2str(bbar)+", $\alpha\in$("+num2str(min(-0.5))+", "+num2str(max(1.5))+")");

plot(eva(1),'rsq',"LineWidth",2);
plot(eva(5),'bsq',"LineWidth",2);
lg=[lg;"Eigenvec 1";"Eigenvec 2"]
legend(lg',"Location","best");


% size_sq23; 
fnts=12; jfm_plt_aid_comm;

subplot(1,4,3)
up=real(evc(:,1)); up=reshape(up.',[4,length(up)/4]);
plot(up,x);

legend("$p$","$u_r$","$u_{\theta}$","$u_z$","Location","best"); grid on; ylabel("z"); title("Eigenvec 1");
fnts=12; jfm_plt_aid_comm;


subplot(1,4,4)
up=real(evc(:,5)); up=reshape(up.',[4,length(up)/4]);
plot(up,x);

legend("$p$","$u_r$","$u_{\theta}$","$u_z$","Location","best"); grid on; ylabel("z"); title("Eigenvec 2");
fnts=12; jfm_plt_aid_comm;

exportgraphics(gcf,"bodespec.eps")

