clc; clear;
% a=importdata("uz.csv");
% a([10,11,47,44,53:55],:)=[];
% a=[-1.1,0;a;1.1,0];

a=importdata("uzE.csv");
a([1,74,39],:)=[];
a=[-1.1,0;a;1.1,0];


clf;
plot(a(:,1),a(:,2));
hold on;
text(a(:,1),a(:,2),num2str([1:length(a)]'))
%%
clf;
x=a(:,1)/1.1; y=a(:,2);
plot(x,y);

hold on;
%
clf;
n=20
xch=-cos(linspace(0,pi,n+1));
T=cos(acos(xch')*[0:n]);
ac=T\interp1(x,y,xch)';
ac(2:2:end)=0;

n2=300;
xch2=(linspace(-1,1,n2));
T2=cos(acos(xch2')*[0:n]);
y2=T2*ac;
plot(xch2,y2)
% stem(xch2,y2)
%%
ac(1:2:end)'