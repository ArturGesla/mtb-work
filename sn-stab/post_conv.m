clc; clear;
%
cd     '/people/gesla/Documents/git/mtb-work/sn-stab'/lorenz_fairgrieve_CN/;
%
ll=importdata("list");
nta=[]; eva=[];
for i=1:length(ll)
    a=load(ll{i}); nta=[nta;a.np ]; eva=[eva;max(a.exp1) ];    i;
end
[a,b]=sort(nta); nta=nta(b); eva=eva(b);
ntaCN=nta; evaCN=eva;

cd     '/people/gesla/Documents/git/mtb-work/sn-stab/lorenz-sn/';
%
ll=importdata("list");
nta=[]; eva=[];
for i=1:length(ll)
    a=load(ll{i}); nta=[nta;a.nt ]; eva=[eva;max(a.ev) ];    i;
end
[a,b]=sort(nta); nta=nta(b); eva=eva(b);
ntaSN=nta; evaSN=eva;

cd /people/gesla/Documents/git/mtb-work/sn-stab/sn-mod-3-lorenz;
%
ll=importdata("list");
nta=[]; eva=[];
for i=1:length(ll)
    a=load(ll{i}); nta=[nta;a.nt ]; eva=[eva;max(a.ev) ];    i;
end
[a,b]=sort(nta); nta=nta(b); eva=eva(b);
ntaCM=nta; evaCM=eva;

cd /people/gesla/Documents/git/mtb-work/sn-stab/lorenz-cheb-coll/;
%
ll=importdata("list");
nta=[]; eva=[];
for i=2:length(ll)
    a=load(ll{i}); nta=[nta;a.nt ]; eva=[eva;max(a.exponents) ];    i;
end
[a,b]=sort(nta); nta=nta(b); eva=eva(b);
ntaCH=nta; evaCH=eva;



%%

plot(nta,eva,'x-');
%%
close all; leg=[]; evEx=evaSN(end);


loglog(ntaCN,abs(evaCN-evEx),'x-'); leg=[leg; "CN "]; hold on; grid on;
loglog(ntaCN,ntaCN.^(-2),'-'); leg=[leg; "2nd order"];
loglog(ntaCN,ntaCN.^(-1),'-'); leg=[leg; "1st order"];

loglog(ntaSN*2,abs(evaSN-evEx),'x-'); leg=[leg; "SN "]; hold on; grid on;
loglog(ntaCM,abs(evaCM-evEx),'x-'); leg=[leg; "CM "]; hold on; grid on;
loglog(ntaCH,abs(evaCH-evEx),'x-'); leg=[leg; "Chebyshev "]; hold on; grid on;


legend(leg);
xlabel("nt*2 or np | 1D DOF"); title("Absolute error in the 0.0465 Fl exponent"); ylabel("abs error");

%%
exportgraphics(gcf,"convergenceStab-2.png")

%% article


close all; leg=[]; evEx=evaSN(end);


loglog(ntaCN,abs(evaCN-evEx),'x-'); leg=[leg; "FD"]; hold on; grid on;
% loglog(ntaCN,ntaCN.^(-1),'-'); leg=[leg; "1st order"];

loglog(ntaSN*2,abs(evaSN-evEx),'o-'); leg=[leg; "Fourier"]; hold on; grid on;
% loglog(ntaCM,abs(evaCM-evEx),'x-'); leg=[leg; "CM "]; hold on; grid on;
loglog(ntaCH,abs(evaCH-evEx),'+-'); leg=[leg; "Chebyshev "]; hold on; grid on;
loglog(ntaCN,20*ntaCN.^(-2),'k-'); leg=[leg; "2nd order"];



legend(leg,"Location","east","Position",[0.654500610783936   0.243501807464159   0.308621007101280   0.227797829115004]);
xlabel("1D dof");% title("Abs err in 0.0465 Fl exponent"); 
ylabel("abslute error");

% jfm_plt_aid;

exportgraphics(gcf,"art3-conv.eps")

%%
close all; clear;
cd     '/people/gesla/Documents/git/mtb-work/sn-stab';
a=load("./lorenz_fairgrieve_CN/stabCNLorenz-21.mat"); plot(real(a.exp1),imag((a.exp1)),'x'); hold on;
a=load("./lorenz-sn/spectrumSNLorenz-10.mat"); plot(real(a.evs),imag((a.evs))/a.om,'o');
a=load("./lorenz-cheb-coll/flnum-20-cheb-coll-lornez.mat"); plot(real(a.exponents),imag((a.exponents)),'+');

xlim([-0.025489432551145   0.065670429159882]); ylim([-10.698402747491688 10.698402747491688]); grid on;
xlabel("$\mu_r$"); ylabel("$\mu_i/ \omega$");
legend("FD","Fourier","Chebyshev","Location","best")

% jfm_plt_aid;


exportgraphics(gcf,"art3-spec.eps")

%%
clc;clear; close all; clrs={[0, 0.4470, 0.7410],[0.8500, 0.3250, 0.0980] 	   ,	[0.9290, 0.6940, 0.1250],[0 0 0]}; mult=1;


a=load("lorenz_fairgrieve_CN/solCNLorenz-31.mat"); iclrs=1;
% iev=3;
% iev=2;
iev=1;
X=reshape(a.u(1:end-1),[3,a.np])'; x=reshape(a.evc(1:end,iev),[3,a.np])'; x1=x; X1=X;
% xp=X+x*mult; plot(xp(:,1),xp(:,2),'.-',"Color",clrs{1}); hold on; a=plot(xp(1,1),xp(1,2),'sq',"Color",clrs{iclrs}); set(a,"HandleVisibility","off"); a=plot(xp(2,1),xp(2,2),'<',"Color",clrs{iclrs}); set(a,"HandleVisibility","off");
xp=X+x*mult; plot3(xp(:,1),xp(:,2),xp(:,3),'.-',"Color",clrs{iclrs}); hold on; a=plot3(xp(1,1),xp(1,2),xp(1,3),'sq',"Color",clrs{iclrs}); set(a,"HandleVisibility","off"); a=plot3(xp(2,1),xp(2,2),xp(2,3),'<',"Color",clrs{iclrs}); set(a,"HandleVisibility","off");

% a=load("lorenz-sn/solSNLorenz-15-67.mat");  iclrs=2;
% a=load("lorenz-sn/solSNLorenz-15-66.mat");  iclrs=2;
a=load("lorenz-sn/solSNLorenz-15-57.mat");  iclrs=2;
 X=a.xBase; x=a.xBasePlusPert-a.xBase; x=x/(x(1,2)./x1(1,2));
% xp=X+x*mult; plot(xp(:,1),xp(:,2),'.-',"Color",clrs{2}); hold on; a=plot(xp(1,1),xp(1,2),'sq',"Color",clrs{iclrs}); set(a,"HandleVisibility","off"); a=plot(xp(2,1),xp(2,2),'<',"Color",clrs{iclrs}); set(a,"HandleVisibility","off");
xp=X+x*mult; plot3(xp(:,1),xp(:,2),xp(:,3),'.-',"Color",clrs{iclrs}); hold on; a=plot3(xp(1,1),xp(1,2),xp(1,3),'sq',"Color",clrs{iclrs}); set(a,"HandleVisibility","off"); a=plot3(xp(2,1),xp(2,2),xp(2,3),'<',"Color",clrs{iclrs}); set(a,"HandleVisibility","off");

% a=load("lorenz-cheb-coll/solChebLorenz-30-3.mat");  iclrs=3;
% a=load("lorenz-cheb-coll/solChebLorenz-30-2.mat");  iclrs=3;
a=load("lorenz-cheb-coll/solChebLorenz-30-1.mat");  iclrs=3;
 X=a.xchBase; x=-a.xchBasePlusPert+a.xchBase; x=x/(x(1,2)./x1(1,2));
% xp=X+x*mult; plot(xp(:,1),xp(:,2),'.-',"Color",clrs{iclrs}); hold on; a=plot(xp(1,1),xp(1,2),'sq',"Color",clrs{iclrs}); set(a,"HandleVisibility","off"); a=plot(xp(2,1),xp(2,2),'<',"Color",clrs{iclrs}); set(a,"HandleVisibility","off");
xp=X+x*mult; plot3(xp(:,1),xp(:,2),xp(:,3),'.-',"Color",clrs{iclrs}); hold on; a=plot3(xp(1,1),xp(1,2),xp(1,3),'sq',"Color",clrs{iclrs}); set(a,"HandleVisibility","off"); a=plot3(xp(2,1),xp(2,2),xp(2,3),'<',"Color",clrs{iclrs}); set(a,"HandleVisibility","off");

%


 iclrs=4;
% xp=X1; plot(xp(:,1),xp(:,2),'k-'); hold on; a=plot(xp(1,1),xp(1,2),'sq',"Color",clrs{iclrs}); set(a,"HandleVisibility","off"); a=plot(xp(2,1),xp(2,2),'<',"Color",clrs{iclrs}); set(a,"HandleVisibility","off");plot(xp(1,1),xp(1,2),'ksq'); plot(xp(2,1),xp(2,2),'k<');
xp=X1; plot3(xp(:,1),xp(:,2),xp(:,3),'.-',"Color",clrs{iclrs}); hold on; a=plot3(xp(1,1),xp(1,2),xp(1,3),'sq',"Color",clrs{iclrs}); set(a,"HandleVisibility","off"); a=plot3(xp(2,1),xp(2,2),xp(2,3),'<',"Color",clrs{iclrs}); set(a,"HandleVisibility","off");

legend("FD","Fourier","Chebyshev","Location","best"); xlabel("x"); ylabel("y");

grid on;    
%
jfm_plt_aid;
%
% exportgraphics(gcf,"art3-lorr24-iev1.eps")
exportgraphics(gcf,"art3-lorr24-iev"+num2str(iev)+"-3d.png","Resolution",150);

