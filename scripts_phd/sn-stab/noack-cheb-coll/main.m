clc; close all; clear; mua=[];
%
% cd     '/people/gesla/Documents/git/mtb-work/sn-stab/sn-mod-noack';
%

% Noack system
np=520; %r=24; 
nt=30; 

% mu=0.04; r=sqrt(mu); gm=1; %gamma
% mu=0.04+it*0.04; r=sqrt(mu); gm=1; %gamma
mu=5; r=sqrt(mu); gm=1; %gamma
neq=3; %np=100; np=np+2;

%init
om1=1;
t=0:2*pi/(np-1):2*pi; u=r*cos(t*om1); v=r*sin(t*om1)/gm; w=r^2*ones(1,length(t));
T=2*pi/om1; x=u; y=v; z=w;

X=[x',y',z'];
%

 X1=X;
tch=cos(0:pi/(np*1-1):pi)'; tch1=tch;

%
tch=(tch+1)/2*t(end); X=interp1(t,X,tch); X=[r*cos(tch*om1),r*sin(tch*om1),r^2*ones(length(tch),1)]; 

y=X; v2=[y;flipud(y(2:end-1,:))]; z=real(fft(v2)./length(v2)); a=z; a(2:end,:)=2*a(2:end,:); %z are cheb coeffs
zcut=z*0; zcut([1:nt,end-nt+2:end],:)=z([1:nt,end-nt+2:end],:); ycut=ifft(zcut).*length(zcut); ycut=ycut(1:end/2+1,:);

%  plot(t,X1,'--'); hold on; set(gca,"ColorOrderIndex",1);
plot(tch,y,':'); hold on; set(gca,"ColorOrderIndex",1);
plot(tch,ycut,'-'); hold on; set(gca,"ColorOrderIndex",1);
fprintf("cheb init approx accuracy: %4.2e\n",norm(y-ycut,"fro"));
grid on; legend("exact x ","exact y ","exact z ","cheb")
%
u=[reshape(real(a(1:nt,:).'),[3*nt,1])];
u(nt*3+1)=2*pi/t(end);
uinit=u;
% semilogy(abs(a(1:2:end,1)))

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
%
% collx=linspace(-1,1,nt)';
collx=-cos(0:pi/(nt*1-1):pi)';
%
close all;
semilogy(abs(reshape(uinit(1:end-1),[3,nt])'),'--'); hold on; grid on;
%
% u=u*0; u(end-3)=1;
for i=1:10
[g,jac]=calculateRhsAndJac(3,nt,u,r,om1,collx);

fprintf("it: %d \t norm(rhs): %4.2e\n",i,norm(g));

if(norm(g)<1e-12)
    break; 
end
u=u-jac\g;
%     u(4:6)'

end
semilogy(abs(reshape(u(1:end-1),[3,nt]))');
% pause;
%
% [U,S,V]=svds(jac(1:end-1,1:end-1),5,'smallest');
% [U,S,V]=svds(jac(1:end,1:end),5,'smallest');
%

% second stability

nt2=nt;
% l3=zeros(3,nt2*2); l3(1,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(1,1)=1; l31=reshape(l3,[3*nt2*2,1])';
% l3=zeros(3,nt2*2); l3(2,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(2,1)=1; l32=reshape(l3,[3*nt2*2,1])';
% l3=zeros(3,nt2*2); l3(3,1:nt2)=ones(1,nt2)*2-mod([0:nt2-1],2)*4; l3(3,1)=1; l33=reshape(l3,[3*nt2*2,1])';

l3=zeros(3,nt2); l3(1,1:nt2)=ones(1,nt2); r31=reshape(l3,[3*nt2,1])';
l3=zeros(3,nt2); l3(2,1:nt2)=ones(1,nt2);  r32=reshape(l3,[3*nt2,1])';
l3=zeros(3,nt2); l3(3,1:nt2)=ones(1,nt2);  r33=reshape(l3,[3*nt2,1])';


jmod=full(jac(1:end-1,1:end-1));%-b*j2hm+b;
ind=nt*3+1-3;
% jmod(ind,:)=l31;
% jmod(ind+1,:)=l32;
% jmod(ind+2,:)=l33;

bmod=jmod*0;
bmod(ind,:)=r31;
bmod(ind+1,:)=r32;
bmod(ind+2,:)=r33;
%
[evc,evs]=eig(full(jmod),full(bmod)); evs=diag(evs); b=abs(evs)<Inf; evs=evs(b); evc=evc(:,b);

evsAnal=[0;(-1-sqrt(1-8*mu))/2;(-1+sqrt(1-8*mu))/2;];
flmult=exp(evsAnal*T);

fprintf("Result || \n");
fprintf("Analytical fl mult:\n");
disp(sort(flmult'));
fprintf("Numerical fl mult:\n");
flnum=1./(1-evs)';
disp(sort(flnum));
fprintf("Diff:\n");
disp(sort(flnum)-sort(flmult'));
% save("flnum-"+num2str(nt)+".mat",'nt','flnum','flmult');
% close all;
% visu
close all;
up=u;
neq=3;
xch=X*0;
for i=0:nt-1    
    xch(:,1)=xch(:,1)+up(i*neq+1).*cos(i*acos(tch1));
    xch(:,2)=xch(:,2)+up(i*neq+2).*cos(i*acos(tch1));
    xch(:,3)=xch(:,3)+up(i*neq+3).*cos(i*acos(tch1));
end

% plot(tch,xch,'-'); hold on; set(gca,"ColorOrderIndex",1); %same as ycut
plot3(xch(:,1),xch(:,2),xch(:,3),'-'); hold on; plot3(xch(1,1),xch(1,2),xch(1,3),'>'); hold on; 
grid on; hold on;

%

% up=u+[V(:,5);0];
up=u+[evc(:,2);0]+[evc(:,1);0];
% up=uinit;
neq=3;
xch=X*0;
for i=0:nt-1    
    xch(:,1)=xch(:,1)+up(i*neq+1).*cos(i*acos(tch1));
    xch(:,2)=xch(:,2)+up(i*neq+2).*cos(i*acos(tch1));
    xch(:,3)=xch(:,3)+up(i*neq+3).*cos(i*acos(tch1));
end

% plot(tch,xch,'-'); hold on; set(gca,"ColorOrderIndex",1); %same as ycut
plot3(xch(:,1),xch(:,2),xch(:,3),'-'); hold on; plot3(xch(1,1),xch(1,2),xch(1,3),'>'); hold on; 
grid on; hold on;

%% post conv
clc; clear; close all;
ll=importdata("list");
fl=[]; nt=[];
for i=1:length(ll)
    a=importdata(ll{i}); a2=a;
    nt=[nt;a.nt]; fl=[fl;sort(a.flnum)];
end

[a,b]=sort(nt); nt=nt(b); fl=fl(b,:);
%%
close all;
semilogy(nt,abs(fl-fl./fl.*sort(a2.flmult)'))
legend(num2str(sort(a2.flmult)))
