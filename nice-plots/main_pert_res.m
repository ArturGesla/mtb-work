clc; clear;
a=read_sunfluidh_data(10); vssH=a.W;  vss=vssH; r=a.yc(1:size(vss,2)); r=r./max(r)*10; z=a.xc(1:size(vss,1)); z=z./max(z); t=a.time;
cc=importdata("resid_L2_Li.d");
%%
 for i=63:63
a=read_sunfluidh_data(i); vss=(abs(vssH-a.W)); aa=a;  r=a.yc(1:size(vss,2)); r=r./max(r)*10; z=a.xc(1:size(vss,1)); z=z./max(z); t=a.time; time=t;
label="log10(u_{\theta}-u_{\theta,steady})";
a=(vss(180,2:end-1)); main_range;
%
vss=log10(abs(vssH-aa.W));
close all;

f=figure(Position=[2200 157 1183 643]); fnts=14;
hold on;
set(f,'defaulttextinterpreter','latex')

% subplot(2,1,1)
tit=("$"+label+" $ $\vert$ $Re_H=1600$ $\vert$ grid 1024x192 $\vert$ $t="+num2str(round(time,2))+"$");
t=tiledlayout(2,3); 
t.TileSpacing = 'tight';
t.Padding = 'loose';
 
nexttile(t,1,[1 3]);
h=pcolor(r,z,vss); h.EdgeAlpha=0;
colorbar(); axis equal; xlim([min(r) max(r)]); ylim([min(z) max(z)]);
xlabel("r"); ylabel("z");  title(tit)
set(gca,"FontSize",fnts,"FontName","Latin Modern Math");
%  subplot(2,1,2)
ax=nexttile(t,4,[1 2]);
h=pcolor(r,z,vss); h.EdgeAlpha=0;
colorbar(); axis equal; xlim([r(i1) r(i2)]); ylim([min(z) max(z)]);
xlabel("r"); ylabel("z"); title(tit)
set(gca,"FontSize",fnts,"FontName","Latin Modern Math");

ax=nexttile(t,6);
semilogy(cc.data(:,1),cc.data(:,2)); title("Residual",FontSize=fnts); hold on;
semilogy([time time],[min(cc.data(2:end,2)) max(cc.data(:,2))],'r')
% xlim([0 1000])
xlabel("t");grid on; grid minor; xlim([1000 3000])
set(gca,"FontSize",fnts,"FontName","Latin Modern Math");

exportgraphics(gcf,"plot"+num2str(i)+".png",'Resolution',300)
 end
%%
f=figure(Position=[2200 202 911 598]); fnts=14;
% hold on;
set(f,'defaulttextinterpreter','latex')

semilogy(c.data(:,1),c.data(:,2)); title("Residual"+"$\vert$ $Re_H=1600$ $\vert$ grid 1024x192",FontSize=fnts)
% xlim([0 1000])
xlabel("t"); ylabel("residual"); grid on; grid minor;
set(gca,"FontSize",fnts,"FontName","Latin Modern Math");
exportgraphics(gcf,'plot.png','Resolution',300)

%%
a=(vss(180,:))
% close all
% plot(a)
% plot(a>max(a)/2)
b=(a>max(a)/2)
c=cumsum(b);
i1=sum(1-cummax(b));
i2=sum(1-cummax(b,"reverse"));
di=270;
ic=(i1+i2)/2;
i1=ic-di/2; i2=ic+di/2;
hold on;
plot(r,a)
plot(r(i1:i2),a(i1:i2)); 
