clc;clear;

r=0:0.1:1; z=r; x=r*0;
[R,Z]=meshgrid(r,z); X=R*0+1;
%%
plot(R,Z,'k'); hold on;
plot(R',Z','k'); hold on;

%%
clc; clf;

z=0.0:0.1:1; th=linspace(0,atan(0.5),20); 
[TH,Z]=meshgrid(th,z);
% X=R*0+1;
%
% X=R; Y=2-R./tan(TH);
% X=2-R.*tan(TH); Y=R;
X=(2-Z).*tan(TH); Y=Z;
plot(X,Y,'k'); hold on;
plot(X',Y','r'); hold on;
%%
mesh(X,Y,Y)

%%


clc;

z=0.0:0.1:1; th=[0.2:0.1:1]*pi/2/2;
[Z,TH]=meshgrid(z,th); 

X=(3-Z).*tan(TH); Y=Z;
plot(X,Y,'k'); hold on;
plot(X',Y','r'); hold on;
