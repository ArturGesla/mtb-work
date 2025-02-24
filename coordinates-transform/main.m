clc;clear;

r=0:0.1:1; z=r; x=r*0;
[R,Z]=meshgrid(r,z); X=R*0+1;
%%
plot(R,Z,'k'); hold on;
plot(R',Z','k'); hold on;

%%
clc;

r=0.5:0.1:1; th=r*pi/2; x=r*0;
[R,TH]=meshgrid(r,th); X=R*0+1;
%
% X=R; Y=2-R./tan(TH);
X=2-R*tan(TH); Y=R;
plot(X,Y,'k'); hold on;
plot(X',Y','r'); hold on;

%%


clc;

z=0.0:0.1:1; th=[0.2:0.1:1]*pi/2/2;
[Z,TH]=meshgrid(z,th); 

X=(3-Z).*tan(TH); Y=Z;
plot(X,Y,'k'); hold on;
plot(X',Y','r'); hold on;
