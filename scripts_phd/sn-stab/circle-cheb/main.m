clc; close all; clear; mua=[];
%
% cd     '/people/gesla/Documents/git/mtb-work/sn-stab/sn-mod-noack';
%

% Lorenz system
np=52; r=24; nt=15; 

T=2*pi;
t=0:2*pi/(np-1):T; t=t';
X=[cos(t),-sin(t),t*0];

 X1=X;
tch=cos(0:pi/(np*1-1):pi)'; tch1=tch;

%
tch=(tch+1)/2*t(end); X=interp1(t,X,tch); % X=[cos(tch),-sin(tch),tch*0];

y=X; v2=[y;flipud(y(2:end-1,:))]; z=real(fft(v2)./length(v2)); a=z; a(2:end,:)=2*a(2:end,:); %z are cheb coeffs
zcut=z*0; zcut([1:nt,end-nt+2:end],:)=z([1:nt,end-nt+2:end],:); ycut=ifft(zcut).*length(zcut); ycut=ycut(1:end/2+1,:);

% plot(t,X1,'--'); hold on; set(gca,"ColorOrderIndex",1);
plot(tch,y,'->'); hold on; set(gca,"ColorOrderIndex",1);
plot(tch,ycut,'-x'); hold on; set(gca,"ColorOrderIndex",1);
% fprintf("cheb init approx accuracy: %4.2e\n",norm(y-ycut,"fro"));
%
u=[reshape(real(a(1:nt,:).'),[3*nt,1])];
u(nt*3+1)=2*pi/T;
% % plot
% neq=3;
% xch=X*0;
% for i=0:nt-1
%     
%     xch(:,1)=xch(:,1)+u(i*neq+1).*cos(i*acos(tch1));
%     xch(:,2)=xch(:,2)+u(i*neq+2).*cos(i*acos(tch1));
%     xch(:,3)=xch(:,3)+u(i*neq+3).*cos(i*acos(tch1));
% end
% 
% plot(tch,xch,'o-'); hold on; set(gca,"ColorOrderIndex",1); %same as ycut
%%
% u=u*0; u(end-3)=1;
for i=1:1
[g,jac]=calculateRhsAndJac(3,nt,u,r);

fprintf("it: %d \t norm(rhs): %4.2e\n",i,norm(g));

if(norm(g)<1e-12)
    break; 
end
u=u-jac\g;
%     u(4:6)'

end
% pause;

%%
% small stab
j2=((jac(1:end-1,1:end-1)));
ev=eigs(j2,1,0.0465)
save("stabSNLorenz-"+num2str(nt)+".mat","ev","nt");

% save("uSN-lorenz-nt9.mat","u")
%%

close all;
clf;

up=u; ntp=nt;
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp);
plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; xpb=xp;
% plot(xp)

%%
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
%%
close all;
clf;

% up=[evc(:,105);0]; ntp=nt; u1=up(1:end-1);
% up=[j2*u2;0]; ntp=nt; 
up=[u3;0]; ntp=nt; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp); xpp=xp;
% xp=xp+xpb;  
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; 

plot(xp)

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