clc; close all; clear; mua=[];
%
% cd     '/people/gesla/Documents/git/mtb-work/sn-stab/sn-mod-noack';
%

% Lorenz system
np=52; r=24;
main_lorenz_ti
%




z=fft(X)./length(X); nt=9; if(mod(nt,2)==0) error("nt even"); end
arr=[1:nt]; %arr(2:2:end)=[];
a1=arr; arr=[arr,length(z)-fliplr(arr(2:end))+2];
zcut=z*0; zcut(arr,:)=z(arr,:);
X2=ifft(zcut)*length(zcut);

close all;
 plot(X); hold on; set(gca,"ColorOrderIndex",1);
plot(X2,':'); hold on
%
z1=zcut(1:nt,:);
om=2*pi/T;
% phase
halfsol = @(phi) real(1i*sum([0:nt-1]'.*z1(:,1).*exp(1i*phi*[0:nt-1]).')); %fully wrong wtf
phi=fzero(halfsol,0);
z1=z1.*exp(i*phi*[0:nt-1]).';
z11=z*0; z11(1:nt,:)=z1; z11(end-nt+2:end,:)=flipud(conj(z1(2:end,:)));
X3=ifft(z11)*length(z11);
set(gca,"ColorOrderIndex",1);plot(X3,'--'); 
%
u=[reshape(real(z1.'),[3*nt,1]);reshape(imag(z1.'),[3*nt,1])];
u(nt*3*2+1)=om/2;
uinit=u;
% load
% load("../lorenz-sn/uSN-lorenz-nt4.mat"); om=u(end); u=reshape(u(1:end-1),[3,4*2]); 
% 
% uinit=reshape(uinit(1:end-1),[3,nt*2]); uinit=uinit*0; uinit(:,1:2:nt)=u(:,1:4); uinit(:,nt+1:2:end)=u(:,5:end);
% u=reshape(uinit,[3*nt*2,1]); u(end+1)=om/2; 
%
%
%  load unt7.mat;
%  load unt13.mat;
%%
for i=1:1
[g,jac]=calculateRhsAndJac(3,nt,u,r);
u=u-jac\g';
fprintf("it: %d \t norm(rhs): %4.2e\n",i,norm(g));

if(norm(g)<1e-12)
    break; 
end
    u(4:6)'

end
%%
svds(jac(1:end-1,1:end-1),5,'smallest')
%
% clf; spy(jac); grid on; grid minor;
% [U,S,V]=svds(jac(1:end,1:end),5,'smallest');s=diag(S);
% [U,S,V]=svds(jac(1:end-1,1:end-1),5,'smallest'); s=diag(S);
% [U,S,V]=svds(jac(1:nt*3,1:nt*3),5,'smallest');
% second stability
%%
nt2=nt;
% l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(1,1)=1; l31=reshape(l3,[3*nt2*2,1])';
% l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(2,1)=1; l32=reshape(l3,[3*nt2*2,1])';
% l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(3,1)=1; l33=reshape(l3,[3*nt2*2,1])';

l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2; l3(1,1)=1; r31=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2; l3(2,1)=1; r32=reshape(l3,[3*nt2*2,1])';
l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2; l3(3,1)=1; r33=reshape(l3,[3*nt2*2,1])';


jmod=full(jac(1:end-1,1:end-1));%-b*j2hm+b;
%  ind=nt*3+1;
%  ind=4:6:(nt-1)*3;ind=[ind,ind+nt*3]
 ind=1*3+1;
% jmod(ind,:)=l31;
% jmod(ind+1,:)=l32;
% jmod(ind+2,:)=l33;

bmod=jmod*0;
bmod(ind,:)=repmat(r31,[1,1]); %r31;
bmod(ind+1,:)=repmat(r32,[1,1]); %r32;
bmod(ind+2,:)=repmat(r33,[1,1]); %r33;

[evc,evs]=eig(full(jmod),full(bmod)); evs=diag(evs); b=abs(evs)<Inf; evs=evs(b); evc=evc(:,b);

fprintf("Result || \n");
fprintf("Numerical fl mult:\n");
disp((1+evs)');
mu=log(evs+1)/pi*u(end);
ev=max(real(mu))
%%
close all;
save("stabCMLorenz-"+num2str(nt)+".mat",'ev','nt');
% load mua.mat; mua=[mua;max(real(mu))]; save('mua.mat','mua');
%% naive stab - Hill
[evc,evs]=eig(full(jac(1:end-1,1:end-1))); evs=diag(evs); %b=abs(evs)<Inf; evs=evs(b); evc=evc(:,b);
%  b=(abs(evs)<Inf); evs=evs(b);
 [a,b]=sort(abs(imag(evs))); 
 fprintf("Numerical fl mult:\n");
disp(exp(2*pi/om*evs(b(1:4))).');
evs=evs(b(1:4)); evc=evc(:,b(1:4)); 
load mua.mat; mua=[mua;max(real(evs))]; save('mua.mat','mua');
%%
% close all;
up=u+[1*evc(:,3);0]; ntp=nt; 
% up=u; ntp=nt; 
% up=[evc(:,2);0]; ntp=nt;  up=up*0; up(4:6)=evc(4:6,2);
% up=[evc(:,3);0]; ntp=nt; 
% up=u; ntp=nt; 
% up=[V(:,5)]; ntp=nt; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';

% up=[V(:,3);0]; ntp=nt; zp=zp+reshape(up(1:(end-1)/2),[3,ntp])'+1*1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';


nz=(np-1)*2-length(zp)*2+1; zp=[zp;zeros(nz,3);conj(flipud(zp(2:end,:)))]; xp=ifft(zp)*(length(zp)); xp=real(xp);

% clf; 
plot3(xp(:,1),xp(:,2),xp(:,3),".-"); grid on; hold on; xpb=xp; plot3(xp(1,1),xp(1,2),xp(1,3),'o');plot3(xp(2,1),xp(2,2),xp(2,3),'>'); plot3(xp(1:end/2,1),xp(1:end/2,2),xp(1:end/2,3),'r');
%  clf; plot(xp); grid on; hold on;
%%
a=load("unt7.mat");
up=a.u; ntp=7; 
zp=reshape(up(1:(end-1)/2),[3,ntp])'+1i*reshape(up((end-1)/2+1:(end-1)),[3,ntp])';
nz=(np-1)*2-length(zp)*2+1; zp=[zp;zeros(nz,3);conj(flipud(zp(2:end,:)))]; xp=ifft(zp)*(length(zp)); xp=real(xp);

%  plot3(xp(:,1),xp(:,2),xp(:,3),'.'); grid on; hold on; xpb=xp; plot3(xp(1,1),xp(1,2),xp(1,3),'o');plot3(xp(2,1),xp(2,2),xp(2,3),'>');
plot(xp)

%%
xp=X; plot3(xp(:,1),xp(:,2),xp(:,3));