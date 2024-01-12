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
open ./lorenz_fairgrieve_CN/main_lor.m
run ./lorenz_fairgrieve_CN/main_lor.m; 
%%
save("dataFD.mat",'X','evc','x','evs');
open ./lorenz-sn/main.m;

run ./lorenz-sn/main.m
save("dataSN.mat",'X4');

%%
load dataFD.mat

%%
plot(X(:,1),X(:,2)); %hold on; plot(X4(:,1),X4(:,2));
%%
plot(X(:,1),X(:,2)); hold on; plot(X(:,1)+x(:,1),X(:,2)+x(:,2)); plot(X(1,1)+x(1,1),X(1,2)+x(1,2),'o');
%%
plot(abs(fft((X(:,1)))),'-x'); hold on; %plot(abs(fft(X4(:,1))),'-o'); 
set(gca,"Yscale","log")
plot(abs(fft((x(1:end-1,1)))),'-x');

%% spec of fd

v=x(1:end-1,:);
z=fft(v)./length(v);
vjr=real(z(1:7,:)); vji=imag(z(1:7,:)); vj=[vjr;vji]; vjj=reshape(vj',[7*3*2,1]);
j=jac(1:end-1,1:end-1);

%% phys space
v2=reshape(vjj,[length(vjj)/2,2]); v2=v2(:,1)+1i*v2(:,2); v21=v2; v2=reshape(v2.',[3,7]).'; v2=[v2;conj(flipud(v2(2:end,:)))];
x2=ifft(v2);%/length(v2);
plot(x(:,1),x(:,2)); hold on; plot(x2(:,1),x2(:,2))

%%
vv=evc(:,39)
v4=vv(1:end/2)+1i*(vv(end/2+1:end))
%%

plot(abs(v21)); hold on; plot(abs(v4)*30); set(gca,"Yscale","log")

abs(v21)./abs(v4)

%%
