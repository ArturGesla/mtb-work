clear; clc; 
om=1;
bbar=0;
R=100;

a=load("../vonkarman-fd-pstag-rotframe/vkrf-np-220.mat");

U=a.u(1:3:end);
V=a.u(2:3:end);
x=a.x;
%%
% D2
z=x;
ii=[]; jj=[]; vv=[];
for i=2:length(x)-1
    zp=z(i); zn=(z(i)+z(i+1))/2; zs=(z(i)+z(i-1))/2;


    ii=[ii;i]; jj=[jj;i]; vv=[vv;(-1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
    ii=[ii;i]; jj=[jj;i]; vv=[vv;(1)/(z(i)-z(i-1))*(1)/(zn-zs)];

    ii=[ii;i]; jj=[jj;i+1]; vv=[vv;(1)/(z(i+1)-z(i))*-(1)/(zn-zs)];

    ii=[ii;i]; jj=[jj;i-1]; vv=[vv;(1)/(z(i)-z(i-1))*-(1)/(zn-zs)];

    D2=-sparse(ii,jj,vv);
end