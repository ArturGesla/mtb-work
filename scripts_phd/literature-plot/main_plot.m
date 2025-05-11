clc;clear; 
gauthiercr=[35000,0.048]; gauthiersr=[62000,0.048];
schouvcr=[35000, 0.048; 23000, 0.07; 18000, 0.085; 16000, 0.1];
schouvsr=[40000, 0.048; 25000, 0.1; 16000, 0.14;];% 16000, 0.1];
poncetcr=[12000, 0.114]; poncetsr=[21000, 0.114];
serresr=[11000, 0.2]; 
rubiosr=[52000, 0.2]; 
gelfgatsr=[22000 0.1; 18000, 0.12; 15000, 0.14; 13000, 0.16; 11500, 0.18; 10000, 0.2];
 daubecr=[300000, 0.1];

tean=[55 128 128]./256;
red=[255 15 0]./256;
%
close all;
f=figure(Position=[2200,200,800,600]); fnts=14;
hold on;
set(f,'defaulttextinterpreter','latex')

plot(serresr(:,1),serresr(:,2),'+',Color=red,LineWidth=2)
plot(rubiosr(:,1),rubiosr(:,2),'+',Color=tean,LineWidth=2)

plot(schouvcr(:,1),schouvcr(:,2),'--',Color=red,LineWidth=2)
plot(schouvsr(:,1),schouvsr(:,2),'--',Color=red,LineWidth=2)

plot(poncetcr(:,1),poncetcr(:,2),'+--',Color=red,LineWidth=2)
plot(poncetsr(:,1),poncetsr(:,2),'+--',Color=red,LineWidth=2)

plot(gauthiercr(:,1),gauthiercr(:,2),'+',Color=tean,LineWidth=2)
plot(gauthiersr(:,1),gauthiersr(:,2),'+',Color=tean,LineWidth=2)

plot(gelfgatsr(:,1),gelfgatsr(:,2),'--',Color=red,LineWidth=2)

% plot(daubecr(:,1),daubecr(:,2),'+',Color=tean,LineWidth=2)

legend("shroud fixed","shroud rotating", "Location",    "best",'Interpreter', 'latex');

% ylim([0.03, 0.22]); xlim([0 350000]); grid on; grid minor;
ylim([0.03, 0.22]); xlim([0 80000]); grid on; grid minor;

% xlabel("$Re_R=\frac{R^2\Omega}{\nu}$", 'Interpreter', 'latex')
xlabel("$Re_R=\frac{R^2\Omega}{\nu}$")
ylabel("$G=\frac{H}{R}$")

% set(gca,'defaulttextinterpreter','latex')

% set(gca, 'FontName', 'CMU Serif');
% set(gca,'DefaultTextFontname', 'CMU Serif')
%    set(gca,'DefaultAxesFontName', 'CMU Serif')


set(gca,"FontSize",fnts,"FontName","Latin Modern Math");

% Create textbox
p=gca().InnerPosition; x0=p(1); y0=p(2); x1=p(3)+x0; y1=p(4)+y0;

xp=60000; yp=0.048;
h=annotation(f,'textbox',...
     [(xp-gca().XLim(1))/diff(gca().XLim)*(x1-x0)+x0 ,(yp-gca().YLim(1))/diff(gca().YLim)*(y1-y0)+y0 , 0.1538, 0.0487],...
    'String',{'SR Gauthier'},'Interpreter','latex',"FitBoxToText","on","FontSize",fnts-2,"LineStyle","none","Rotation",   0);

xp=28000; yp=0.035;
h=annotation(f,'textbox',...
     [(xp-gca().XLim(1))/diff(gca().XLim)*(x1-x0)+x0 ,(yp-gca().YLim(1))/diff(gca().YLim)*(y1-y0)+y0 , 0.1538, 0.0487],...
    'String',{'CR Gauthier'},'Interpreter','latex',"FitBoxToText","on","FontSize",fnts-2,"LineStyle","none","Rotation",   0);

xp=26000; yp=0.1;
h=annotation(f,'textbox',...
     [(xp-gca().XLim(1))/diff(gca().XLim)*(x1-x0)+x0 ,(yp-gca().YLim(1))/diff(gca().YLim)*(y1-y0)+y0 , 0.1538, 0.0487],...
    'String',{'SR Schouveiler'},'Interpreter','latex',"FitBoxToText","on","FontSize",fnts-2,"LineStyle","none","Rotation",   -45);

% xp=13000; yp=0.08;
% h=annotation(f,'textbox',...
%      [(xp-gca().XLim(1))/diff(gca().XLim)*(x1-x0)+x0 ,(yp-gca().YLim(1))/diff(gca().YLim)*(y1-y0)+y0 , 0.1538, 0.0487],...
%     'String',{'CR Schouveiler'},'Interpreter','latex',"FitBoxToText","on","FontSize",fnts-2,"LineStyle","none","Rotation",   -45);
% 
% xp=9000; yp=0.1;
% h=annotation(f,'textbox',...
%      [(xp-gca().XLim(1))/diff(gca().XLim)*(x1-x0)+x0 ,(yp-gca().YLim(1))/diff(gca().YLim)*(y1-y0)+y0 , 0.1538, 0.0487],...
%     'String',{'CR \& SR Poncet'},'Interpreter','latex',"FitBoxToText","on","FontSize",fnts-2,"LineStyle",":","Rotation",   0,"BackgroundColor",[1,1,1]);

xp=12000; yp=0.18;
h=annotation(f,'textbox',...
     [(xp-gca().XLim(1))/diff(gca().XLim)*(x1-x0)+x0 ,(yp-gca().YLim(1))/diff(gca().YLim)*(y1-y0)+y0 , 0.1538, 0.0487],...
    'String',{'SR Gelfgat'},'Interpreter','latex',"FitBoxToText","on","FontSize",fnts-2,"LineStyle","none","Rotation",   -70);

xp=11000; yp=0.2;
h=annotation(f,'textbox',...
     [(xp-gca().XLim(1))/diff(gca().XLim)*(x1-x0)+x0 ,(yp-gca().YLim(1))/diff(gca().YLim)*(y1-y0)+y0 , 0.1538, 0.0487],...
    'String',{'SR Serre'},'Interpreter','latex',"FitBoxToText","on","FontSize",fnts-2,"LineStyle","none","Rotation",   0);

xp=51000; yp=0.2;
h=annotation(f,'textbox',...
     [(xp-gca().XLim(1))/diff(gca().XLim)*(x1-x0)+x0 ,(yp-gca().YLim(1))/diff(gca().YLim)*(y1-y0)+y0 , 0.1538, 0.0487],...
    'String',{'SR Lopez'},'Interpreter','latex',"FitBoxToText","on","FontSize",fnts-2,"LineStyle","none","Rotation",   0);

% xp=300000; yp=0.1;
% h=annotation(f,'textbox',...
%      [(xp-gca().XLim(1))/diff(gca().XLim)*(x1-x0)+x0 ,(yp-gca().YLim(1))/diff(gca().YLim)*(y1-y0)+y0 , 0.1538, 0.0487],...
%     'String',{'CR Daube'},'Interpreter','latex',"FitBoxToText","on","FontSize",fnts-2,"LineStyle","none","Rotation",   0);
% 
exportgraphics(gcf,'plot.png','Resolution',300)
