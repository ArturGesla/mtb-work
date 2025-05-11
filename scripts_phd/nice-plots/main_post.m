clc; clear;
a=read_sunfluidh_data(10); vss=a.W; r=a.yc; z=a.xc; t=a.time;
c=importdata("resid_L2_Li.d");
%%
for i=120:130
% i=0;
a=read_sunfluidh_data(i); v=a.W-vss; r=a.yc; z=a.xc; t=a.time;
%
f=figure();
mesh(r,z,log(abs(v)+v)); set ( gca, 'ydir', 'reverse' ); xlabel("r"); ylabel("z"); title("Vth| Reh=1600 | t = " +num2str(t)+" | mesh 1024x192 nu");
% colorbar(); colormap(gca,parula(16))
axes('Position',[.6 .6 .3 .3])
box on
semilogy(c.data(:,1),c.data(:,2)); hold on; semilogy([t,t],[1e-20,1e-0]); grid on;
% saveas(gca,"vth"+num2str(i)+".jpg")
exportgraphics(gcf,"vth"+num2str(i)+".jpg",'Resolution',400)
close(f);
end
%%

semilogy(c.data(:,1),c.data(:,2)); title("Residual")
hold on;
semilogy(c.data(:,1),c.data(:,3)); xlabel("t"); ylabel("residual"); grid on; grid minor;
legend("L2 norm","Linf norm")
%%

semilogy(c.data(:,1),c.data(:,2)); title("Residual")
hold on;
% semilogy(c.data(:,1),c.data(:,3)); 
xlabel("t"); ylabel("residual"); grid on; grid minor;
legend("L2 norm","Linf norm")
