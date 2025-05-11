mu1=[-0.25:0.001:0.25];
mu3=[-0.25:0.001:0];
mu5=[0:0.01:0.25];
mu7=[-0.28:0.01:0.0];
delta=1+4*mu1;
r1=sqrt((1+sqrt(delta))/2);
delta=1+4*mu3;
r3=sqrt((1-sqrt(delta))/2);

close all;
f=figure(); fnts=14; set(f,"Visible","on"); 
set(f,"Position",[ 850         594        300 120]);
set(gca,"Position",[0.1300    0.200    0.7750    0.8150])
% f=figure(Position=[2200 202 300 300]); fnts=14;
% hold on;
set(f,'defaulttextinterpreter','latex')
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
%

hold on;
al=1;
% plot(mu1,r1*al,'k')
% plot(mu3,r3*al,'k--')
% plot(mu5,mu5*0,'k--')
% plot(mu7,mu7*0,'k')



mu2=[-0.28:0.001:0.25];
mu4=[0:0.01:0.25];

r2=real(sqrt(mu2));
plot(mu2,r2*al,'-k');
plot(mu4,mu4*0,'--k');


ylim([-.01 1.2*al]);
xlim([-0.3 0.3]);
% axis off;
set(gca,"XTick",[0])
% set(gca,"XTick",[])
% set(gca,"XColor","none")
% set(gca,"X
set(gca,"XTickLabel","Re_c")
set(gca,"Color",'none')
set(gca,'TickLength',[0 0])
% set()
ax=gca;
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';
text(0,-.01,"$Re_c$",'HorizontalAlignment', 'center','VerticalAlignment','top','FontSize',14,'FontName','Helvetica')
% ax.XAxis.Label.Visible = 'on';

% set(gca,"FontSize",fnts,"FontName","Latin Modern Math");
set(gca,"FontSize",10,"FontName","Helvetica");
% set(gca,"LabelFontSizeMultiplier",5)
% title("lol")

% annotation('textarrow',[0.2 0.2],[0 1],'String',' Growth ','FontSize',13,'Linewidth',1)
% plot([-0.28 -0.28],[0 1.1],'-k')
% text(-0.33,0.8,"A",'FontSize',14,'FontName','Helvetica')

plot([-0.28 -0.28],[0 1.1],'-k')
text(-0.33,0.8,"A",'FontSize',14,'FontName','Helvetica')

exportgraphics(gcf,"plotsup.eps","BackgroundColor","none")
% print('plotsub','-dsvg','-r300')
% print('plotsub','-deps','-r300')
%%

close all;
f=figure("Position",[670   677   344   287]);
set(f,"OuterPosition",[670   677   344   372])

set(f,'defaulttextinterpreter','latex')
hold on;
leg=[];


a=importdata("top.dat"); a=importdata("200498.grapdir/control-600-160-2300-nu.dat");
plot(a(:,1),a(:,2),'-k'); leg=[leg;"chaotic rolls"];
a=importdata("bot.dat"); a=importdata("200499.grapdir/control-600-160-2300-nu.dat");
plot(a(:,1),a(:,2),'-k'); leg=[leg;"chaotic rolls"];
a=importdata("latest.out");
plot(a(:,1),a(:,2),'-r'); leg=[leg;"chaotic rolls"];
% plot(bottomBranch(1,:),bottomBranch(2,:),'--or'); leg=[leg;"edge"];


% grid on; grid minor; 
xlim([0 1000]); 
% xlim([1600 3400]); 
ylim([0 8.4]);
% title("Bifurcation diagram rotor-stator H/R=0.1");
set(gca,"FontSize",12,"FontName","Latin Modern Math");
set(gca,"Position",[0.1300    0.200    0.7750    0.8150])

xlabel("time"); ylabel("$\vert \vert \omega_{pert} \vert \vert _{L2} $");
% legend(leg,"Location","northwest");
% exportgraphics(gcf,"bifurcationDiagramPub.png","Resolution",300)
set(gca,"Color",[255,253,250]./256)
set(gca,"Position",[0.1300    0.200    0.7750    0.8150])


exportgraphics(f,"aaplotTE.eps","BackgroundColor",[255,253,250]./256)
print("bifurcationDiagramPub.png",'-dpng','-r300')
% set(gcf,"Color",[255,253,250]./256)
% print('plotTE','-depsc','-r300')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');

%% subplot

close all;
f=figure("Position",[670   677   704   287]);
% set(f,"OuterPosition",[670   677   344   372])

set(f,'defaulttextinterpreter','latex')
leg=[];


subplot(1,2,1)
hold on;

plot([2925.47],[0],"ksq"); leg=[leg;"Re_c 600-160"];% text([2925.47],[0.1],"Re_c")
% errorbar(bottomBranch(1,:),bottomBranch(2,:),bottomBranch(5,:),'--r'); leg=[leg;"std"];

% errorbar(topBranch(1,:),topBranch(2,:),(topBranch(3,:)-topBranch(4,:))/2,'-ok'); leg=[leg;"chaotic rolls"];
% errorbar(bottomBranch(1,:),bottomBranch(2,:),(bottomBranch(3,:)-bottomBranch(4,:))/2,'--ob'); leg=[leg;"edge"];

plot(topBranch(1,:),topBranch(2,:),'-ok'); leg=[leg;"chaotic rolls"];
plot(bottomBranch(1,:),bottomBranch(2,:),'--or'); leg=[leg;"edge"];


% grid on; grid minor; 
xlim([1600 3400]); ylim([0 8.4]);

set(gca,"FontSize",12,"FontName","Latin Modern Math");

xlabel("$Re=\Omega H^2/\nu$"); ylabel("$\vert \vert \omega_{pert} \vert \vert _{L2} $");
% legend(leg,"Location","northwest");
% exportgraphics(gcf,"bifurcationDiagramPub.png","Resolution",300)
set(gca,"Color",[255,253,250]./256)


subplot(1,2,2)
hold on;

a=importdata("top.dat"); a=importdata("200498.grapdir/control-600-160-2300-nu.dat");
plot(a(:,1),a(:,2),'-k'); leg=[leg;"chaotic rolls"];
a=importdata("bot.dat"); a=importdata("200499.grapdir/control-600-160-2300-nu.dat");
plot(a(:,1),a(:,2),'-k'); leg=[leg;"chaotic rolls"];
a=importdata("latest.out");
plot(a(:,1),a(:,2),'-r'); leg=[leg;"chaotic rolls"];
% plot(bottomBranch(1,:),bottomBranch(2,:),'--or'); leg=[leg;"edge"];



% grid on; grid minor; 
xlim([0 800]); 
% xlim([1600 3400]); 
ylim([0 8.4]);
% title("Bifurcation diagram rotor-stator H/R=0.1");
set(gca,"FontSize",12,"FontName","Latin Modern Math");
xlabel("time");
set(gca,"Color",[255,253,250]./256)

exportgraphics(f,"aaplotBoth.eps","BackgroundColor",[255,253,250]./256)
% print("bifurcationDiagramPub.png",'-dpng','-r300')
% set(gcf,"Color",[255,253,250]./256)
% print('plotTE','-depsc','-r300')
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');