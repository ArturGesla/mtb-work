clc; close all; clear;
%
np=52; r=24; nt=15;
main_lorenz_ti
%
z=fft(X); 
arr=[1:nt]; a1=arr; arr=[arr,length(z)-fliplr(arr(1:end-1))+1];
zcut=z*0; zcut(arr,:)=z(arr,:);
X2=ifft(zcut);
%
% close all;
% plot(X); hold on;
% plot(X2); hold on;
%
z1=zcut(a1,:)./np;
om=2*pi/T;
an=angle(sum(z1(2,1))); %this could be better
% z1=z1.*exp(-an*1i*[0:nt-1]');
angle(sum(z1(:,1)));
% zcut(arr,:)=zcut(arr,:).*exp(an*1i*[arr-1]');
% X3=ifft(zcut);
% plot(X3); hold on;
% x=real(z1(1,1)+z1(2,1)*exp(i*om*t)+z1(3,1)*exp(2*i*om*t)+z1(2,1)'*exp(-i*om*t)+z1(3,1)'*exp(-2*i*om*t));
x=real(z1(1,1)+z1(2,1)*exp(1i*om*t)+z1(2,1)'*exp(-1i*om*t));
plot(x,'-x')
%
% z1=zcut(a1,:);
u=[reshape(real(z1.'),[3*nt,1]);reshape(imag(z1.'),[3*nt,1])];
u(nt*3*2+1)=om;
u0=u;
%%
close all;
%
% r=r+1e-6
for i=1:16
[g,jac]=calculateRhsAndJac(3,nt,u,r);
u=u-jac\g';
fprintf("it: %d \t norm(rhs): %4.2e \t gend: %4.2e\n",i,norm(g),g(end));
if(norm(g)<1e-10)
    break;
end
end
%%
% semilogy(abs(u(1:3:nt*3)+u(nt*3+1:3:end-1)),'x-'); hold on;
% semilogy(abs(u(2:3:nt*3)+u(nt*3+2:3:end-1)),'x-'); hold on;
% semilogy(abs(u(3:3:nt*3)+u(nt*3+3:3:end-1)),'x-'); hold on;
% grid on;
% %%
% semilogy(abs(u(1:3:nt*3)),'xb-'); hold on;
% semilogy(abs(+u(nt*3+1:3:end-1)),'xb--'); hold on;
% 
% semilogy(abs(u(2:3:nt*3)+u(nt*3+2:3:end-1)),'x-'); hold on;
% semilogy(abs(u(3:3:nt*3)+u(nt*3+3:3:end-1)),'x-'); hold on;
% grid on;
% %% small stab
% j2=((jac(1:end-1,1:end-1)));
% ev=eigs(j2,1,0.0465)
% save("stabSNLorenz-"+num2str(nt)+".mat","ev","nt");
% 
% % save("uSN-lorenz-nt9.mat","u")
% %%
% 
% close all;
% clf;
% 
% up=u; ntp=nt;
% zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
% xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xpb=xp;
% % plot(xp)

%
close all;
% evsM={}; leg=[];
rank(full(jac(1:end-1,1:end-1)))
j2=(full(jac(1:end-1,1:end-1)));
    B=j2*0; 
for i=[1:nt*3,nt*3+4:nt*3*2]
    B(i,i)=1; 
end
% 
[evc,evs]=eig(j2,B); evs=diag(evs);
% evsM{end+1}=evs; leg=[leg;"nt = "+num2str(nt)];

% om=2*pi/T;
om=u(end);
 plot(real(evs),imag(evs)/om,'+'); grid on; hold on;
 text(real(evs),imag(evs)/om,num2str([1:length(evs)]'));

 b=(abs(evs)<Inf); evs=evs(b);
 [a,b]=sort(abs(imag(evs))); 
 fprintf("Numerical fl mult:\n");
disp(exp(2*pi/om*evs(b(1:3))).');
save("spectrumSNLorenz-"+num2str(nt)+".mat","evs","nt","om");

%%
close all;
clf;

% up=[evc(:,105);0]; ntp=nt; u1=up(1:end-1);
% up=[j2*u2;0]; ntp=nt; 
up=[u0]; ntp=nt; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp); xpp=xp;
% xp=xp+xpb;  
plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xBase=xp;

%%
iev=57; up=u+[evc(:,iev);0]; ntp=nt; u1=up(1:end-1);
% up=[j2*u2;0]; ntp=nt; 
% up=[u]; ntp=nt; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp); xpp=xp;
% xp=xp+xpb;  
plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xBasePlusPert=xp;

% plot(xp)

% x for  cheb
xp=[xp;xp(1,:)];
t=linspace(0,2*pi/u(end),length(xp));
% save("xforcheb"+num2str(r)+".mat","xp","t",'r','u','nt');

%
save("solSNLorenz-"+num2str(nt)+"-"+num2str(iev)+".mat","xBasePlusPert","xBase","u",'r');

%% expo counterpart
ll=length(xp); sigma=evs(105)*1;
t=0:T/ll:T-T/ll;
xp=xpp.*(exp(sigma*t)');

z2=fft(xp)./length(xp);
u2=z2(1:nt,:); 
u2=[reshape(real(u2.'),[3*nt,1]);reshape(imag(u2.'),[3*nt,1])]; 
%%
clf; ip=2; 
% plot(u2(ip:3:end)); hold on;
u3=j2*u2; plot(u3(ip:3:end))