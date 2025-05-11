clear; clc;
syms r z Re;
ur = 3*r^3*z^2;
uz=-4*r^2*z^3;
p=r*0;


fr=-1/Re*(+24*r*z^2+6*r^3)+3*r^5*z^4;
fz=1/Re*(+24*r^2*z+16*z^3)+24*r^4*z^5;


divu=diff(ur,r)+ur/r+diff(uz,z);

% p=1/Re*(12*r^2*z^2+6/4*r^4)-27/6*r^6*z^4+24/6*r^6*z^4;

% momr=ur*diff(ur,r)+uz*diff(ur,z)+diff(p,r)-1/Re*(1/r*diff(r*diff(ur,r),r)-ur/r/r+diff(diff(ur,z),z));
momr=ur*diff(ur,r)+uz*diff(ur,z)+diff(p,r)-1/Re*(1/r*diff(r*diff(ur,r),r)-ur/r/r+diff(diff(ur,z),z))-fr
momz=ur*diff(uz,r)+uz*diff(uz,z)-1/Re*(1/r*diff(r*diff(uz,r),r)+diff(diff(uz,z),z))-fz
%%

cd /home/gesla/Documents/git/rotst2/build;
addpath     '/home/gesla/Documents/git/rotst2/scripts/source_for_mtb';

a=importdata("u20-100-120.dat");
a=loadrs2(a,100,120);
%%
mesh(a.xc,a.zc,a.v)

%%
mesh(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1),'EdgeColor','r');

hold on;
mesh(a.xu,a.zc(1:end-1,1:end-1),3*a.xu.^3.*a.zc(1:end-1,1:end-1).^2,'EdgeColor','b');
%%

clf; hold on; grid on;
errArr=[]; nxArr=[];

nx=4; nz=5;
system("rotst -mode 3 -nx "+num2str(nx)+" -nz "+num2str(nz)+" -nnewt  10 -lx 1 -lz 1 -regEps 201 0");
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
mean(a.zc((nz+3)/2,:))
plot(a.xu(1,:),a.u((nz+3)/2,1:end-1),'-x')
errArr=[errArr;a.xu(1,:).^3*0.5^2*3-a.u((nz+3)/2,1:end-1)];
nxArr=[nxArr; nx*nz];
%

nx=nx*2; nz=nz*2-1;
system("rotst -mode 3 -nx "+num2str(nx)+" -nz "+num2str(nz)+" -nnewt  10 -lx 1 -lz 1 -regEps 201 0");
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
mean(a.zc((nz+3)/2,:))
plot(a.xu(1,:),a.u((nz+3)/2,1:end-1),'-o')
errArr=[errArr;a.xu(1,1:2:end).^3*0.5^2*3-a.u((nz+3)/2,1:2:end-1)];
nxArr=[nxArr; nx*nz];


nx=nx*2; nz=nz*2-1;
system("rotst -mode 3 -nx "+num2str(nx)+" -nz "+num2str(nz)+" -nnewt  10 -lx 1 -lz 1 -regEps 201 0");
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
mean(a.zc((nz+3)/2,:))
plot(a.xu(1,:),a.u((nz+3)/2,1:end-1),'-+')
errArr=[errArr;a.xu(1,1:4:end).^3*0.5^2*3-a.u((nz+3)/2,1:4:end-1)];
nxArr=[nxArr; nx*nz];

nx=nx*2; nz=nz*2-1;
system("rotst -mode 3 -nx "+num2str(nx)+" -nz "+num2str(nz)+" -nnewt  10 -lx 1 -lz 1 -regEps 201 0");
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
mean(a.zc((nz+3)/2,:))
plot(a.xu(1,:),a.u((nz+3)/2,1:end-1),'-+')
errArr=[errArr;a.xu(1,1:8:end).^3*0.5^2*3-a.u((nz+3)/2,1:8:end-1)];
nxArr=[nxArr; nx*nz];

nx=nx*2; nz=nz*2-1;
system("rotst -mode 3 -nx "+num2str(nx)+" -nz "+num2str(nz)+" -nnewt  10 -lx 1 -lz 1 -regEps 201 0");
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
mean(a.zc((nz+3)/2,:))
plot(a.xu(1,:),a.u((nz+3)/2,1:end-1),'-+')
errArr=[errArr;a.xu(1,1:16:end).^3*0.5^2*3-a.u((nz+3)/2,1:16:end-1)];
nxArr=[nxArr; nx*nz];

%%
loglog(sqrt(nxArr),errArr(:,2:end-1),'-x')
hold on;
loglog(sqrt(nxArr),sqrt(nxArr).^(-2))
%%
(errArr(1:end-2,:)-errArr(2:end-1,:))./(errArr(2:end-1,:)-errArr(3:end,:))
%%
clc; clf; 
errArr=[]; nxArr=[];


nx=4; nz=5;
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
% mesh(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1));
du=a.u(1:end-1,1:end-1)-3*a.xu.^3.*a.zc(1:end-1,1:end-1).^2;
mesh(a.xu,a.zc(1:end-1,1:end-1),du);
errArr=[errArr;max(max(abs(du)))];
nxArr=[nxArr; nx];


nx=nx*2; nz=nz*2-1;
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
% mesh(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1));
du=a.u(1:end-1,1:end-1)-3*a.xu.^3.*a.zc(1:end-1,1:end-1).^2;
mesh(a.xu,a.zc(1:end-1,1:end-1),du);
errArr=[errArr;max(max(abs(du)))];
nxArr=[nxArr; nx];


nx=nx*2; nz=nz*2-1;
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
% mesh(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1));
du=a.u(1:end-1,1:end-1)-3*a.xu.^3.*a.zc(1:end-1,1:end-1).^2;
mesh(a.xu,a.zc(1:end-1,1:end-1),du);
errArr=[errArr;max(max(abs(du)))];
nxArr=[nxArr; nx];


nx=nx*2; nz=nz*2-1;
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
% mesh(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1));
du=a.u(1:end-1,1:end-1)-3*a.xu.^3.*a.zc(1:end-1,1:end-1).^2;
mesh(a.xu,a.zc(1:end-1,1:end-1),du);
errArr=[errArr;max(max(abs(du)))];
nxArr=[nxArr; nx];


nx=nx*2; nz=nz*2-1;
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
% mesh(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1));
du=a.u(1:end-1,1:end-1)-3*a.xu.^3.*a.zc(1:end-1,1:end-1).^2;
mesh(a.xu,a.zc(1:end-1,1:end-1),du);
errArr=[errArr;max(max(abs(du)))];
nxArr=[nxArr; nx];
%%

nx=24; nz=29;
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
% mesh(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1));
du=a.u(1:end-1,1:end-1)-3*a.xu.^3.*a.zc(1:end-1,1:end-1).^2;
mesh(a.xu,a.zc(1:end-1,1:end-1),du);
errArr=[errArr;max(max(abs(du)))];
nxArr=[nxArr; nx];


nx=48; nz=57;
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
% mesh(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1));
du=a.u(1:end-1,1:end-1)-3*a.xu.^3.*a.zc(1:end-1,1:end-1).^2;
mesh(a.xu,a.zc(1:end-1,1:end-1),du);
errArr=[errArr;max(max(abs(du)))];
nxArr=[nxArr; nx];
%%
clf;
loglog(nxArr,errArr); hold on;
loglog(nxArr,nxArr.^(-1))

%%


nx=4; nz=5;
% nx=nx*2; nz=nz*2-1;
% nx=nx*2; nz=nz*2-1;
% nx=nx*2; nz=nz*2-1;
% nx=nx*2; nz=nz*2-1;

a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);
% mesh(a.xu,a.zc(1:end-1,1:end-1),a.u(1:end-1,1:end-1));
contour(a.xu(2:end,:),a.zc(2:end-1,1:end-1),a.u(2:end-1,1:end-1),[0:0.1:1].^3*3,'-k');
% du=a.u(1:end-1,1:end-1)-3*a.xu.^3.*a.zc(1:end-1,1:end-1).^2;
% mesh(a.xu,a.zc(1:end-1,1:end-1),du);
% nxArr=[nxArr; nx];

%%
nx=40; nz=40;
a=importdata("u20-"+num2str(nx)+"-"+num2str(nz)+".dat");
a=loadrs2(a,nx,nz);



%%
arr=[
%     20, 4.73332474424090821e-01;
    40,4.74073179047988347e-01;
    80,4.74272609395770650e-01;
    160,4.74324166312154072e-01];


arr=[
%     20, 4.73332474424090821e-01;
    40,4.97932064344595582e-01;
    80,4.86167458083175463e-01;
    160,4.80262610742940632e-01];

arr=[
%     20, 4.73332474424090821e-01;
    40,4.74813708911960841e-01;
    80,4.74457861250872359e-01;
%         4.74370486415251702e-01
    160,4.74370486415251702e-01];


arr=[
%     20, 4.73332474424090821e-01;
%     40,6.17755232969368961e-01;
    80,6.17350176594676947e-01;
%         4.74370486415251702e-01
    160,6.17247768602027014e-01;
    320,6.17222015038383698e-01];



% exact=sqrt(9/8/5);
exact=sqrt(16/6/7);
or=(arr(1,2)-arr(2,2))./(arr(2,2)-arr(3,2));
loglog(arr(:,1),abs(arr(:,2)-exact),'-x'); hold on; grid on;
loglog(arr(:,1),arr(:,1).^(-2),'-');
loglog(arr(:,1),arr(:,1).^(-1),'-');
