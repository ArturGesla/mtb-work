%%
clf;
% close all;
%
hold on;
n=4; x=-1:0.001:1;
% t=0:0.001:pi; x=cos(t);
acos(x);
y=acos(x);
z=cos(n*y);
% plot(z,'-x')
plot(x,z,'-'); 
grid on;