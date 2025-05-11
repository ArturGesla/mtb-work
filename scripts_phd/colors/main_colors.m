arr=[];
c='63ccc5'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr=[arr;a];
c='176885'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr=[arr;a];
c='001c2e'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr=[arr;a];

plt=["#001C30","#176B87","#64CCC5","#DAFFFB"]; arr5=[];
% c='001C30'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
% c='176B87'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
% c='64CCC5'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
% c='DAFFFB'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];


% c='112D4E'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
% c='3F72AF'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
% c='DBE2EF'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
% c='F9F7F7'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];

%not bad
c='1B262C'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
c='0F4C75'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
c='3282B8'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
c='BBE1FA'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];

% c='2C3333'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
% c='395B64'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
% c='A5C9CA'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];
% c='E7F6F2'; a=[string(c(1:2)),c(3:4),c(5:6)]; a=[hex2dec(a(1)),hex2dec(a(2)),hex2dec(a(3))]./256; arr5=[arr5;a];

% plot([1,2],[1,2],"Color",a);

% plot3(arr(:,1),arr(:,2),arr(:,3))
n=3;
arr2=[linspace(arr(1,1),arr(3,1),n*2)',linspace(arr(1,2),arr(3,2),n*2)',linspace(arr(1,3),arr(3,3),n*2)'];
arr3=[linspace(arr(1,1),arr(2,1),n+1),linspace(arr(2,1),arr(3,1),n+1);
    linspace(arr(1,2),arr(2,2),n+1),linspace(arr(2,2),arr(3,2),n+1);
    linspace(arr(1,3),arr(2,3),n+1),linspace(arr(2,3),arr(3,3),n+1)]';
arr3(n+1:n+2,:)=[];
arr3=flipud(arr3)+0.2;
%%
close all;
x=-10:0.1/2:10;
[x,y]=meshgrid(x,x); colormap(arr3)

% figure()
% pcolor(x,y,x); %shading interp; 
% colormap(arr);
% figure()
% pcolor(x,y,sin(10*atan2(y,x)).*cos(sqrt(x.^2+y.^2)));
% colormap(arr5);shading interp; 

figure()
% pcolor(x,y,cos(sqrt(x.^2+y.^2)));
pcolor(x,y,sin(10*atan2(y,x)).*cos(sqrt(x.^2+y.^2)));
colormap(arr3);shading interp; 


%%
cd C:\Users\Artur\Documents\work
a=importdata("u1800.dat.base"); nx=602;nz=162;neq=4; a=a(3:neq:end); v=reshape(a,[nz,nx]); vb=v;
% a=importdata("u1800.dat"); nx=602;nz=162;neq=4; a=a(3:neq:end); v=reshape(a,[nz,nx]); v=v-vb;
%%
% close all;
% figure()
% pcolor(x,y,sin(10*atan2(y,x)).*cos(sqrt(x.^2+y.^2)));
% colormap(arr2);shading interp; 

figure()
% pcolor(x,y,cos(sqrt(x.^2+y.^2)));
pcolor(v(2:end-1,2:end-1));
colormap(arr3);shading interp; colorbar();
