addpath("/people/gesla/Documents/git/rotst2/scripts/source_for_mtb")
%%
x=0:1/10:1-1/10; y=exp(0*x).*(sin(1*x*2*pi)+0.2*sin(2*x*2*pi)+0.1*sin(3*x*2*pi));
z=fft(y)/length(y);
% z=fft(y);
close all; plot([0:length(x)-1],abs(z(1:end)),'x-'); 
% hold on; plot(imag(z),'x-');


y=exp(1*x).*(sin(1*x*2*pi)+0.2*sin(2*x*2*pi)+0.1*sin(3*x*2*pi));
z=fft(y)/length(y);
% z=fft(y);
hold on; plot([0:length(x)-1],abs(z(1:end)),'x-'); 

%%
close all;
plot(real(z(2:end)))

z2=[z(1:end/2+1),zeros(1,40),z(end/2+1:end)]
x2=length(z2)*ifft(z2);

plot([0:1/(length(x2)-1):1],x2,'x-');
hold on;
plot(x,y,'-o')

%%
% open ./lorenz_fairgrieve_CN/main_lor.m
run ./lorenz_fairgrieve_CN/main_lor.m; 
%
save("dataFD.mat");
%open ./lorenz-sn/main.m;

run ./lorenz-sn/main.m
save("dataSN.mat");

%%
clear;
fd=load("dataFD.mat");
sn=load("dataSN.mat");

%% base flow

vv=sn.u(1:end-1); vjj=vv;
v2=reshape(vjj,[length(vjj)/2,2]); v2=v2(:,1)+1i*v2(:,2); v21=v2; v2=reshape(v2.',[3,7]).'; v2=[v2;conj(flipud(v2(2:end,:)))];
z2=v2;
X2=real(ifft(z2)*length(z2));
fd.X=fd.X(1:end-1,:);


plot(fd.X(:,1),fd.X(:,2),'-x'); hold on;% text(fd.X(:,1),fd.X(:,2),num2str([1:length(fd.X(:,2))]')); 
plot(fd.X(1,1),fd.X(1,2),'ro');plot(fd.X(2,1),fd.X(2,2),'r>'); %hold on; plot(X4(:,1),X4(:,2));

% dx=(fd.X(:,1)-X2(1,1)'); dy=(fd.X(:,2)-X2(1,2)'); dz=(fd.X(:,3)-X2(1,3)'); d=dx.^2+dz.^2+dz.^2; [a,b]=min(d);
% fdX=fd.X;
% fdX=fdX([b:length(fdX),1:b-1],:);
% plot(fdX(1,1),fdX(1,2),'ro');plot(fdX(2,1),fdX(2,2),'r>'); %hold on; plot(X4(:,1),X4(:,2));

plot(X2(:,1),X2(:,2));plot(X2(1,1),X2(1,2),'bo');plot(X2(2,1),X2(2,2),'b>'); %hold on; plot(X4(:,1),X4(:,2));
% %%
% close all;
% plot(X(:,1),X(:,2)); hold on; plot(X(:,1)+x(:,1),X(:,2)+x(:,2)); plot(X(1,1)+x(1,1),X(1,2)+x(1,2),'o');
% %%
% plot(abs(fft((X(:,1)))),'-x'); hold on; %plot(abs(fft(X4(:,1))),'-o'); 
% set(gca,"Yscale","log")
% plot(abs(fft((x(1:end-1,1)))),'-x');
% 
%% spec of fd

v=fd.x(1:end-1,:);
z=fft(v)./length(v);
z=z(1:7,:);
z=[z;conj(flipud(z(2:end,:)))];
xx=ifft(z)*length(z);
plot(0:1/length(xx):1-1/length(xx),xx); hold on; plot(0:1/length(v):1-1/length(v),v); % w punktach probkowania dokladny overlap
%% phys space
% v2=reshape(vjj,[length(vjj)/2,2]); v2=v2(:,1)+1i*v2(:,2); v21=v2; v2=reshape(v2.',[3,7]).'; v2=[v2;conj(flipud(v2(2:end,:)))];
% x2=ifft(v2);%/length(v2);
% plot(x(:,1),x(:,2)); hold on; plot(x2(:,1),x2(:,2))

%% spec of sn
j=sn.jac(1:end-1,1:end-1);
[evc,evs]=eig(full(j)); evs=diag(evs);
vv=evc(:,39); vjj=vv; %36+om 37-om 32 33
v2=reshape(vjj,[length(vjj)/2,2]); v2=v2(:,1)+1i*v2(:,2); v21=v2; v2=reshape(v2.',[3,7]).'; v2=[v2;conj(flipud(v2(2:end,:)))];
z2=v2;

vv=evc(:,38); vjj=vv; %38 0.0465
v2=reshape(vjj,[length(vjj)/2,2]); v2=v2(:,1)+1i*v2(:,2); v21=v2; v2=reshape(v2.',[3,7]).'; v2=[v2;conj(flipud(v2(2:end,:)))];
z38=v2;

vv=evc(:,29); vjj=vv; %29 -13
v2=reshape(vjj,[length(vjj)/2,2]); v2=v2(:,1)+1i*v2(:,2); v21=v2; v2=reshape(v2.',[3,7]).'; v2=[v2;conj(flipud(v2(2:end,:)))];
z29=v2;


%% terrible ideas ahead
% vv=evc(:,37)+evc(:,36); vjj=vv; %36 37 32 33
% v2=reshape(vjj,[length(vjj)/2,2]); v2=v2(:,1)+1i*v2(:,2); v21=v2; v2=reshape(v2.',[3,7]).'; v2=[v2;conj(flipud(v2(2:end,:)))];
% za=v2;
% 
% j=full(j);
% j(end,:)=0; j(end,end/2+1:end)=reshape([[0:6]',zeros(7,2)]',[1 21]);
% b=eye(42); b(end,end)=0;
% [evc2,evs2]=eig(full(j),b); evs2=diag(evs2);
% 
% vv=evc2(:,38); vjj=vv; %36 37 32 33
% v2=reshape(vjj,[length(vjj)/2,2]); v2=v2(:,1)+1i*v2(:,2); v21=v2; v2=reshape(v2.',[3,7]).'; v2=[v2;conj(flipud(v2(2:end,:)))];
% zbc=v2;
%% porownanie spec

%% abs
clf;
plot(abs(z./norm(z)),'-x'); 
hold on;
plot(abs(z2./norm(z2)),'-+'); %evc normalised wrt the const
plot(abs(z29./norm(z29)),'-x'); 
plot(abs(z38./norm(z38)),'-x'); 
set(gca,"Yscale","log")
% plot(abs(za./norm(za)),'-x'); 
% z2=za;

%%
% z=z./norm(z); z2=z2./norm(z2);
xx=ifft(z)*length(z);
% xx2=(ifft(z2.*[exp(4i*(0:12)/13*2*pi)].')*length(z2)); 
xx2=(ifft(z2)*length(z2)); disp(norm(imag(xx2))); xx2=real(xx2);
xx29=(ifft(z29)*length(z29)); disp(norm(imag(xx29))); xx29=real(xx29);
xx38=(ifft(z38)*length(z38)); disp(norm(imag(xx38))); xx38=real(xx38);
% plot(0:1/length(xx):1-1/length(xx),xx); hold on; %plot(0:1/length(v):1-1/length(v),v); % w punktach probkowania dokladny overlap
% plot(0:1/length(xx2):1-1/length(xx2),xx2);
% plot(1-1/length(xx2):-1/length(xx2):0,-xx2);
% legend(num2str([1:6]'))

%
clf
plot(xx(:,1),xx(:,2)); hold on; plot(v(:,1),v(:,2)); 
hold on;
%  xx2=-xx2*1.5/1.2; %evs up to a asclar 
plot(xx2(:,1),xx2(:,2))
% xx29=-xx29;
plot(xx29(:,1),xx29(:,2))
plot(xx38(:,1),xx38(:,2))
% plot(-xx2(:,1),-xx2(:,2))

%%
clf;

xp=X2; plot3(xp(:,1),xp(:,2),xp(:,3)); hold on;
% xp=X2+xx2; plot3(xp(:,1),xp(:,2),xp(:,3));
% xp=X2+xx38; plot3(xp(:,1),xp(:,2),xp(:,3));
xp=X2+xx29; plot3(xp(:,1),xp(:,2),xp(:,3));

T=2*pi/sn.u(end); t=0:T/length(X2):T-T/length(X2);
% xp=X2+xx38.*exp(20*evs(38)*t)'; plot3(xp(:,1),xp(:,2),xp(:,3));
xp=X2+xx29.*exp(evs(29)*t)'; plot3(xp(:,1),xp(:,2),xp(:,3));
% view([0,0,1])
%% ft comparison
xx29=xx38;
evs(29)=evs(38);
z29=fft(xx29);
xe29=xx29.*exp(evs(29)*t)';
ze29=fft(xe29);

plot(abs(z29)); hold on; plot(abs(ze29));
% plot([xx29,xe29]);
%% mappings

zv29=reshape(z29(1:7,:).',[21,1]); zv29=[real(zv29);imag(zv29)];
a=zv29; b=j*zv29;
%%
zve29=reshape(ze29(1:7,:).',[21,1]); zve29=[real(zve29);imag(zve29)];
a=zve29; b=j*zve29;

