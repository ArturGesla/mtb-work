%%
clc;
clf;
clear; 
load lorenz-sn/spectrumSNLorenz-28.mat;
% up=[evc(:,105);0]; ntp=nt; u1=up(1:end-1);
% up=[j2*u2;0]; ntp=nt; 
up=[u]; ntp=nt; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;zeros(1000,3);conj(flipud(zp(2:end,:)))]);
% *((length(zp)-1)*2+1);
xp=xp*length(xp);
xp=real(xp); xpp=xp;
% xp=xp+xpb;  
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xBase=xp;

XF=xp;
UF=zp;
T=2*pi/u(end);
mu2=evs2(2);
%

up=[evc(:,2)]; ntp=nt; 
zp=reshape(up(1:(end)/2),[3,ntp])'+1i*reshape(up((end)/2+1:(end)),[3,ntp])';
xp=ifft([zp;zeros(1000,3);conj(flipud(zp(2:end,:)))]);
% *((length(zp)-1)*2+1);
xp=xp*length(xp);
xp=real(xp); xpp=xp;


XFprime=xp;
UFprime=zp;
%%

load lorenz-cheb-coll/flnum-150-160-cheb-coll-lornez.mat;
UC=reshape(u(1:end-1),[3 nt])';
UCprime=reshape(evc(1:end,2),[3 nt])';
%%
% clf; set(gcf,"Position",[     510   248   925   522]);
clf; set(gcf,"Position",[      510   155   925   615]);
addpath ../../rotst2/scripts/source_for_mtb/;

tiledlayout(2,3);

nexttile;
plot([0:length(XF)-1]./(length(XF)-1),XF);
xlabel('t'); 
legend('x','y','z',"Location","northwestoutside"); grid on;title("Base solution");
fnts=12; jfm_plt_aid_comm; ylim([-75 225])

nexttile;
semilogy([0:length(UF)-1],abs(UF),'x'); grid on; ylim([1e-6 1e3]); xlim([0 27]); title("mag. Fourier coefficients");
 jfm_plt_aid_comm;

nexttile;
semilogy([0:length(UC)-1],abs(UC),'.'); grid on; ylim([1e-9 1e3]); xlim([0 150]); title("mag. Chebyshev coefficients");
jfm_plt_aid_comm;

%pert
nexttile;
t=linspace(0,T,1+length(XFprime)); t=t(1:end-1);
plot([0:length(XFprime)-1]./(length(XFprime)-1),XFprime.*exp(t*mu2)');
xlabel('t'); 
legend("x'","y'","z'","Location","northwestoutside"); grid on;title("Least stable eigenvector");
 jfm_plt_aid_comm; 
%

nexttile;
semilogy([0:length(UFprime)-1],abs(UFprime),'x'); grid on;
ylim([1e-4 1e1]);
xlim([0 27]); title("mag. Fourier coefficients");
 jfm_plt_aid_comm;

nexttile;
semilogy([0:length(UCprime)-1],abs(UCprime),'.'); grid on; 
% ylim([1e-9 1e3]); 
xlim([0 150]); title("mag. Chebyshev coefficients");
 jfm_plt_aid_comm;

%%
exportgraphics(gcf,"spectralContent.pdf")


