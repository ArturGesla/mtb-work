clc; clear;
%%
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
for i=1:length(ll)
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

