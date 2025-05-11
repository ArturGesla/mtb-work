clc; clear;
a=read_sunfluidh_data(10); vss=a.W;  r=a.yc(1:size(vss,2)); r=r./max(r)*10; z=a.xc(1:size(vss,1)); z=z./max(z); t=a.time;
a=read_sunfluidh_data(120); vss=vss-a.W;  r=a.yc(1:size(vss,2)); r=r./max(r)*10; z=a.xc(1:size(vss,1)); z=z./max(z); t=a.time; time=t;
label="u_{\theta}";

%%
close all;

f=figure(Position=[2200 202 911 598]); fnts=14;
hold on;
set(f,'defaulttextinterpreter','latex')



% subplot(2,1,1)
tit=("$"+label+" $ $\vert$ $Re_H=1600$ $\vert$ grid 1024x192 $\vert$ $t="+num2str(round(time,2))+"$");
t=tiledlayout(2,1); 
t.TileSpacing = 'tight';
t.Padding = 'loose';
 
nexttile;
h=pcolor(r,z,vss); h.EdgeAlpha=0;
colorbar(); axis equal; xlim([min(r) max(r)]); ylim([min(z) max(z)]);
xlabel("r"); ylabel("z");  title(tit)
set(gca,"FontSize",fnts,"FontName","Latin Modern Math");
%  subplot(2,1,2)
ax=nexttile;
h=pcolor(r,z,vss); h.EdgeAlpha=0;
colorbar(); axis equal; xlim([1.5 4.5]); ylim([min(z) max(z)]);
xlabel("r"); ylabel("z"); title(tit)
set(gca,"FontSize",fnts,"FontName","Latin Modern Math");

exportgraphics(gcf,'plot.png','Resolution',300)

%%
f=figure(Position=[2200 202 911 598]); fnts=14;
% hold on;
set(f,'defaulttextinterpreter','latex')

semilogy(c.data(:,1),c.data(:,2)); title("Residual"+"$\vert$ $Re_H=1600$ $\vert$ grid 1024x192",FontSize=fnts)
% xlim([0 1000])
xlabel("t"); ylabel("residual"); grid on; grid minor;
set(gca,"FontSize",fnts,"FontName","Latin Modern Math");
xlabel("t"); ylabel("$\vert\vert \Phi^n-\Phi^{n-1}\vert\vert$ ")
exportgraphics(gcf,'plot.png','Resolution',300)
