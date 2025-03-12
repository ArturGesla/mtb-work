clc; clear;
a=importdata("opf.csv"); %hub wall save


%
Y=(reshape(a.data(:,2),[65,129]));
Z=(reshape(a.data(:,3),[65,129]));
U=(reshape(a.data(:,4),[65,129]));
V=(reshape(a.data(:,5),[65,129]));
W=(reshape(a.data(:,6),[65,129]));
P=(reshape(a.data(:,7),[65,129]));
%
Y=Y(:,3:end);
Z=Z(:,3:end);
W=W(:,3:end);
% mesh(Z,Y,W)
mesh(Z,Y,W)

%%
% system("rotst -mode 3 -nx 64 -nz 128 -lx 1 -lz 4 -shroudRot 0 -contipath 1 10 -regEps 100 0")
% system("rotst -mode 3 -nx 8 -nz 32 -lx 1 -lz 4 -shroudRot 0 -contipath 1 100 -regEps 100 0")
system("rotst -mode 3 -nx 16 -nz 64 -lx 1 -lz 4 -shroudRot 0 -contipath 1 100 -regEps 100 0")
%%
% a=importdata("u10-64-128.dat"); 
a=importdata("u100-8-32.dat"); 
addpath     '/home/gesla/Documents/git/rotst2/scripts/source_for_mtb';
% a=loadrs2(a,64,128);
a=loadrs2(a,8,32);
%%
mesh(a.xc,a.zc,a.w)
%%

cd ../pipeFull/
%%
clc; 
clear;
a=importdata("opf.csv"); %hub wall save
%

Y=(reshape(a.data(1:18*33,2),[18,33]));
Z=(reshape(a.data(1:18*33,3),[18,33]));
W=(reshape(a.data(1:18*33,6),[18,33]));

mesh(Z,Y,W)
hold on;
 clf;
a=importdata("u100-8-32.dat"); 
a=loadrs2(a,8,32);
mesh(a.zc,a.xc,a.w)
%
hold on;
vi=interp2(Z(1,:),Y(:,1),W,a.zc(2:end-1,1),a.xc(1,2:end-1))';
mesh(a.zc(2:end-1,1),a.xc(1,2:end-1),interp2(Z(1,:),Y(:,1),W,a.zc(2:end-1,1),a.xc(1,2:end-1)))

dv=a.w(2:end-1,2:end-1)-vi;
%%
mesh(dv)

%%

clc; 
clear;
a=importdata("opf32.csv"); %hub wall save
%

X=(reshape(a.data(1:35*65,1),[35,65]));
Z=(reshape(a.data(1:35*65,3),[35,65]));
W=(reshape(a.data(1:35*65,6),[35,65]));
mesh(Z,X,W)
% mesh(Y)

hold on;
%  clf;
a=importdata("u100-16-64.dat"); 
a=loadrs2(a,16,64);
mesh(a.zc,a.xc,a.w)
%%
% hold on;
vi=interp2(Z(1,:),X(:,1),W,a.zc(2:end-1,1),a.xc(1,2:end-1))';
mesh(a.zc(2:end-1,1),a.xc(1,2:end-1),interp2(Z(1,:),X(:,1),W,a.zc(2:end-1,1),a.xc(1,2:end-1)))

dv=a.w(2:end-1,2:end-1)-vi;

%%

mesh(dv(:,1:10))
