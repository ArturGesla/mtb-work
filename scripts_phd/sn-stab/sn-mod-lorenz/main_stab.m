%% new method - remap to half mode
nt2=(nt-1)*2+1;
u2=[]; u2(nt2*3*2+1)=om/2;

z=reshape(u(1:(end-1)/2),[3,nt])'+1i*reshape(u((end-1)/2+1:(end-1)),[3,nt])';
z2=reshape(u2(1:(end-1)/2),[3,nt2])'+1i*reshape(u2((end-1)/2+1:(end-1)),[3,nt2])';

z2(1:2:end,:)=z;

u2(1:end-1)=[real(reshape(z2.',[nt2*3,1])); imag(reshape(z2.',[nt2*3,1]))]; u2=u2.';
% u2 should solve g

[g,jac2]=calculateRhsAndJac(3,nt2,u2, mu,gm); norm(g)
% [g,jac2]=calculateRhsAndJac(3,nt2,u2); norm(g)

j2hm=jac2(1:end-1,1:end-1);
[evc2,evs2]=eig(full(j2hm)); evs2=diag(evs2);
[evc,evs]=eig(full(jac(1:end-1,1:end-1))); evs=diag(evs);
%
clf;
plot(evs,'x'); grid on;
hold on; 
% plot(evs2,'+'); grid on;
% text(real(evs2),imag(evs2),num2str([1:length(evs2)]')); 

ll3=zeros(1,length(evs)); ind1=3; ll3(1:3+ind1)=1; ll3(nt*3+1:nt*3+4+ind1)=1;
b=diag(ll3);

[evc,evs]=eig(full(jac(1:end-1,1:end-1)),b); evs=diag(evs);
plot(real(evs),imag(evs),'o'); grid on;
xlim([-1 0.2])

evsAnal=[0;(-1-sqrt(1-8*mu))/2;(-1+sqrt(1-8*mu))/2;];
% plot(real(evsAnal),imag(evsAnal),'>k');
%%
ll=mod(floor((0:nt2*3-1)./3),2);
ll2=[ll,ll];
% b=zeros(length(evs2)); b(4:6,4:6)=eye(3);
b=diag(ll2);
%%
%% 
% [evc3,evs3]=eig(full(j2hm),b); evs3=diag(evs3); plot(evs3,'sq'); grid on;
% bb=abs(evs3)<Inf; evc3=evc3(:,bb); evs3=evs3(bb); 
% [evc3,evs3]=eig(full(j2hm-b*j2hm+b),b); evs3=diag(evs3); plot(evs3,'sq'); grid on;



%
% spy(j2hm-b*j2hm+b);
% spy(b)
% 
% grid on; set(gca,"XTick",[1:3:nt2*3*2]);set(gca,"YTick",[1:3:nt2*3*2])

%
l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(1,1)=1; l31=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(2,1)=1; l32=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(3,1)=1; l33=reshape(l3,[3*nt2*2,1])';

l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2; l3(1,1)=1; r31=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2; l3(2,1)=1; r32=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2; l3(3,1)=1; r33=reshape(l3,[3*nt2*2,1])';

%
jmod=j2hm;%-b*j2hm+b;
jmod(28,:)=l31;
jmod(29,:)=l32;
jmod(30,:)=l33;

bmod=jmod*0;
bmod(28,:)=r31;
bmod(29,:)=r32;
bmod(30,:)=r33;

%
[evc,evs]=eig(full(jmod),full(bmod)); evs=diag(evs); b=abs(evs)<Inf; evs=evs(b); evc=evc(:,b);

%%
close all;
clf;

up=[evc(:,1);0]; ntp=nt2; u1=up(1:end-1);
% up=[j2*u2;0]; ntp=nt; 
% up=[u2;0]; ntp=nt; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
xp=ifft([zp;conj(flipud(zp(2:end,:)))])*((length(zp)-1)*2+1);xp=real(xp); xpp=xp;
% xp=xp+xpb;  
% plot3(xp(:,1),xp(:,2),xp(:,3)); grid on; hold on; 

plot(xp)
