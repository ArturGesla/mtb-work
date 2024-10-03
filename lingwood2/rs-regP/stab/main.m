 clc; clear;
 
 addpath C:\Users\Artur\Documents\GitHub\rotst2\scripts\source_for_mtb;
%%
% a=load("../rs-np-132-k-0.313-L-32.mat");
% a=load("../rs-np-200-k-0.313-L-31.6228.mat");
% a=load("../rs-np-162-k-0.313-reh-1000.mat");
% a=load("../rs-np-162-k-0.313-reh-600.mat");
% a=load("../rs-np-162-k-0.313-reh-500.mat");
addpath ..
a=get_rs(1000,160);

data=a;
x=a.x;
u=a.u*0;
U=a.u;
z=x;
zw=[2*z(1)-z(2),z(1:end-1),2*z(end-1)-z(end-2)];
zc=(zw(1:end-1)+zw(2:end))/2;
zw=zw(2:end);
k=a.k;
% omega=0+0.04i;
% % beta=bbar=beta/R;
% bbar=0.126;
% R=1;
% alpha=1; 

% [g,jac0,jac1,jac2]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);

% R=515; bbar=0.0117;
% beta=0; bbar=beta/R;

up=reshape(a.u,[4,length(a.u)/4])';
plot(up(:,1:end),a.x,'k-'); 

% aa=load('data2D.mat'); aa=aa.data;
% up=reshape(aa.u,[4,length(aa.u)/4])';
% hold on; plot(up(:,1:end),aa.x,'-'); 

%%

evaM=[];
%%
 
R=27.4;  bbar=0.1152; beta=bbar*R; 
omegaa=[-0.1:0.01/4/4:0.1]+0.01i/2/2/4/4*0; alpha=0;

R=40; bbar=-0.0316; shift=0.1i; 
R=  40.127834854531294;
omegaa=[-0.05:0.01/10:0.05]+1i*0.001/4; alpha=0;

R=50; bbar=0;
omegaa=[-0.05:0.01/10:0.05]+1i*1e-3; alpha=0;

% R=R/sqrt(k); bbar=beta/R;  omegaa=omegaa.*k^(3/2); %ar=ar.*sqrt(k);


%% temporal
eva=[]; im=[];
% for i=1:120;%120
    for i=1:length(omegaa);%120
omega=omegaa(i);

[g,jac0,jac1,jac2,jacom]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
% tic; [evc,evs]=polyeig(jac0,jac1,jac2); ev=(evs); toc; 
tic; [evc,evs]=polyeigs2(jac0+jacom*omega,jac1,jac2,20,0.2); ev=diag((evs)); toc; 
eva=[eva;ev]; im=[im;ev*0+omega];
% omega=omega+0.01/2/2/4%/2
disp(omega)
i
end
% evaM=[evaM,eva];
%%
save("eva-omline-k-"+num2str(k)+"-l-"+num2str(length(x))+".mat",'eva','k','x');
%%
clf; grid on; 
a=load("eva-omline-k-1.mat"); plot(a.eva,'r.'); hold on;
a=load("eva-omline-k-0.313.mat"); plot(a.eva/k^(1/2),'bsq');
%% alplot

clf;
plot(eva,'k.'); hold on;
% plot(eva,'k.'); hold on; text(real(eva),imag(eva),num2str(im))
% plot(evaM,'.'); hold on;
% plot(eva,'x'); hold on;
%  plot(ev,'x'); hold on; text(real(ev),imag(ev),num2str([1:length(ev)]'))
% plot(ev,'o');
% xlim([0 0.35]);
ylim([-0.2 0.2]*2); 
grid on;
% xlim([-0.1 1.5]); ylim([-1 1]); grid on;
% xlim([-0.1 0.5]*10); ylim([-0.5 0.5]*10); grid on;
% xlim([-0.1 0.5]*40); ylim([-0.5 0.5]*40)

xlim([-1 1]); ylim([-0.6 0.6 ]); xlabel("$a_r$"); ylabel("$a_i$"); 
title("Re="+num2str(R)+", $\beta$="+num2str(bbar)+", $\omega\in$("+num2str(min(omegaa))+", "+num2str(max(omegaa))+")");
size_sq23; fnts=10; jfm_plt_aid_comm;
% exportgraphics(gcf,"p4-lw97fig6-1.eps")%lw97 fig 6

% pbaspect([5 1 1]); xlabel("a_r"); ylabel("a_i"); 
% % fnts=12; jfm_plt_aid_comm;
% % exportgraphics(gcf,"spatialbranches-c.eps")
% % exportgraphics(gcf,"spatialbranches-b.eps")
% exportgraphics(gcf,"spatialbranches-a.eps")

% save("eva-"+num2str(length(a.x))+".mat",'eva');

%%
clf;
iev=5;
up=reshape(real(evc(:,iev)),[4,length(x)])';
plot(up(:,[1:3]),zc,'x-'); grid on;
hold on;
plot(up(:,4),zw,'x-');
% plot(up(:,1:1),x,'x-'); grid on;
title("ev "+num2str(iev)+":"+num2str(ev(iev),'%4.2e'))
iev=iev+1
ylim([0 max(zw)])
%%
evaM=[];
%% spatial
omega=0;
shift=0.1i+0.1;
shift=0;
ai=0;
eva=[];
z=[];

% ai=-0.5:0.1:0.2;
ar=-0.05:0.01:0.45; 


% ai=-0.2:0.01:0;
% ar=0.1:1e-2:0.3;

% R=560; bbar=0.10;
% % ai=[0];  ar=-0.05:0.01:0.45; 
% ai=[-0.2:0.01:0];  ar=0:0.01:0.4; 

R=27.4; bbar=0.1152; beta=bbar*R;
ai=[0];  ar=-0.5:0.01:1.5; %ar=0.5;
% ai=[-0.5:0.01:0.5];  ar=0:0.01:1; 
R=R/sqrt(k); bbar=beta/R;ar=ar*sqrt(k); shift=shift/k^(3/2); %ar=ar./sqrt(k);

R=40; bbar=-0.0316; shift=0.1i; 
ai=[0];  ar=-0.5:0.01:1.5; %ar=0.5;
ar=0:0.01/2:1; ai=0:0.01:0.1;
% ar=0.1812; ai=0.0259;
R=  40.127834854531294; ai=0; %ai=0.025936806686725;


R=50; bbar=0; shift=0.1+0.1i;
% ai=[0]+0.00;  ar=-0.5:0.01:1.5; ar(51)=[];%ar=0.5;
% ar=-0.3:0.01:0.3; ai=-0.1:0.01:0.1; ar(31)=[];
% ar=0.2:0.01:0.4; ai=0:0.01:0.15; 
% ar=0.2:0.01:0.3; ai=0.06;
ar=0.01:0.01/2:0.4; ai=-0.05:0.01/4/2:0.05; %reh=500

% R=50; bbar=0; ar=0.01:0.01:1; ai=0; %reh=600
%
zRe=[];
% for ire=10:10:300
%     R=ire;
for ii=1:length(ai) 
for ir=1:length(ar) 
    alpha=ar(ir)+1i*ai(ii);
     [g,jac0,jac1,jac2,jacom]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
    jac=jac0+alpha*jac1+alpha.^2*jac2;
%     [ev]=eigs(jac,-jacom,30,shift); eva=[eva,ev];
    [evc,ev]=eigs(jac,-jacom,30,shift); ev=diag(ev); eva=[eva;ev];
%     disp(alpha);
fprintf("ir %d\t ii %d\n",ir,ii);
    [~,b]=max(imag(ev));
    z(ir,ii)=ev(b);
    end
end
% zRe=[zRe, z];
% disp(ire);
% end
% evaM=[evaM,eva];
%%
clf;
plot(eva,'k.'); grid on; hold on;
plot(ev,'r.'); hold on; text(real(ev),imag(ev),string(num2str([1:length(ev)]')));

% plot(evaM,'.'); grid on; hold on;

% %%
% hold on; axis equal;
% plot(shift,'rx')
% 
% xlim([-0.08 0.04]); ylim([-0.08 0.01 ]); xlabel("$\omega_r$"); ylabel("$\omega_i$"); 
% title("Re="+num2str(R)+", $\beta$="+num2str(bbar)+", $\alpha\in$("+num2str(min(ar))+", "+num2str(max(ar))+")");
% size_sq23; fnts=10; jfm_plt_aid_comm;
% exportgraphics(gcf,"p4-lw97fig6-2.eps")%lw97 fig 6
%%
% save("eva-areal-k-"+num2str(k)+"-l-"+num2str(length(x))+"+.mat",'eva','k','x');
save("evars-areal-k-"+num2str(k)+"-np-"+num2str(length(x))+"-L-"+num2str(x(end))+".mat",'eva','k','x','evc','ev');
%%
clf; grid on; hold on;
% a=load("eva-areal-k-1-np-130-L-30.mat"); plot(a.eva,'r.'); hold on;
a=load("../../bode-regP/stab/eva-areal-k-0.313-np-130-L-30.mat"); plot(a.eva/a.k^(3/2),'r+'); 
a=load("evars-areal-k-0.313-np-132-L-32.mat"); plot(a.eva/a.k^(3/2),'bo'); 
%%
clf;
% a=load("eva-areal-k-1.mat"); plot(a.eva,'k.'); hold on;
a=load("eva-areal-k-1-np-140-L-30.mat"); plot(a.eva,'rx'); hold on;
a=load("eva-areal-k-1-np-130-L-30.mat"); plot(a.eva,'b+'); hold on;
a=load("eva-areal-k-1-np-130-L-60.mat"); plot(a.eva,'gsq'); hold on;

%% evc comp
clf;
% a=load("eva-areal-k-1-np-100-L-30.mat"); plot(a.eva,'b+'); hold on;  ev=a.eva; %text(real(ev),imag(ev),string(num2str([1:length(ev)]')));
% a=load("eva-areal-k-1-np-200-L-60.mat"); plot(a.eva,'gsq'); hold on; ev=a.eva; %text(real(ev),imag(ev),string(num2str([1:length(ev)]')));
a=load("eva-areal-k-1-np-1600-L-480.mat"); plot(a.eva,'rx'); hold on; ev=a.eva; text(real(ev),imag(ev),string(num2str([1:length(ev)]')));

% a=load("eva-areal-k-1-np-100-L-30.mat"); up=reshape(abs(a.evc(:,9)),[4,length(a.x)])'; plot(up(:,[1:4]),a.x,'r-'); grid on; hold on;
% a=load("eva-areal-k-1-np-200-L-60.mat"); up=reshape(abs(a.evc(:,18)),[4,length(a.x)])'; plot(up(:,[1:4]),a.x,'b-'); grid on;
%% z cont
clf;
contour(ar,ai,imag(z)',40); title("Imaginary part of most unstable $\omega(\alpha)$.");
hold on; contour(ar,ai,imag(z)',[0 0],'-k');
hold on; contour(ar,ai,imag(z)',[1 1]*7.8e-4,'k--');

hold on; contour(ar,ai,imag(z)',[1 1]*5.98e-4,'k:');

% contour(ar,ai,real(z)',400); title("Real part of most unstable $\omega(\alpha)$.");

% surf(ar,ai,imag(z)');
xlabel("$\alpha_r$");
ylabel("$\alpha_i$");
fnts=12; jfm_plt_aid_comm; 
size_sq; colorbar(); colormap(hsv(8));

%%

exportgraphics(gcf,"saddlepoint-rs-ss.pdf");


%% precise saddle

clc;
% x=[0.25;-0.05];
% x=[0.2;0.2];
% x=[0.2;0.02];
 bbar=0; R=50; x=[0.25;    0.05];% reh=200
  bbar=0; R=50; x=[0.15;    -0.01];% reh=200

%  x=[0.173;    -0.00135];% reh=500
%  x=[0.244;    8.56e-3];% reh=500
xa=[x];
eps=1e-6;
zh=@(x,y) imag(imagOmega(x,y,bbar,R,data));

%%

d2fdx2=(zh(x(1)+eps,x(2))-2*zh(x(1),x(2))+zh(x(1)-eps,x(2)))/eps^2;
d2fdxdy=(zh(x(1)+eps,x(2)+eps)+zh(x(1)-eps,x(2)-eps)-zh(x(1)+eps,x(2)-eps)-zh(x(1)-eps,x(2)+eps))/eps^2/4;
d2fdy2=(zh(x(1),x(2)+eps)-2*zh(x(1),x(2))+zh(x(1),x(2)-eps))/eps^2;

jac=[d2fdx2,d2fdxdy; d2fdxdy, d2fdy2];
g=[(zh(x(1)+eps,x(2))-zh(x(1)-eps,x(2)))/2/eps;(zh(x(1),x(2)+eps)-zh(x(1),x(2)-eps))/2/eps];
x=x-jac\g;
xa=[xa,x];
det(jac);
norm(g);
% imagOmega(x(1),x(2),bbar,R,data)

fprintf("ar: %4.2e\tai: %4.2e\tR: %4.2f\tnorm(g): %4.2e\tomi: %4.2e\n",x(1),x(2),R,norm(g),zh(x(1),x(2)));

%% 3dof Re

clc; x=[0.25;    0.05; 50];
% x=[0.244;    8.56e-3; 50];% reh=500
% x=[0.173;    -0.00135; 50];% reh=500

xa=[x];
eps=1e-6;

zh=@(x,y,R) imag(imagOmega(x,y,bbar,R,data));

%%


d2fdx2=(zh(x(1)+eps,x(2),x(3))-2*zh(x(1),x(2),x(3))+zh(x(1)-eps,x(2),x(3)))/eps^2;
d2fdxdy=(zh(x(1)+eps,x(2)+eps,x(3))+zh(x(1)-eps,x(2)-eps,x(3))-zh(x(1)+eps,x(2)-eps,x(3))-zh(x(1)-eps,x(2)+eps,x(3)))/eps^2/4;
d2fdy2=(zh(x(1),x(2)+eps,x(3))-2*zh(x(1),x(2),x(3))+zh(x(1),x(2)-eps,x(3)))/eps^2;

d2fdxdb=(zh(x(1)+eps,x(2),x(3)+eps)+zh(x(1)-eps,x(2),x(3)-eps)-zh(x(1)+eps,x(2),x(3)-eps)-zh(x(1)-eps,x(2),x(3)+eps))/eps^2/4;
d2fdydb=(zh(x(1),x(2)+eps,x(3)+eps)+zh(x(1),x(2)-eps,x(3)-eps)-zh(x(1),x(2)+eps,x(3)-eps)-zh(x(1),x(2)-eps,x(3)+eps))/eps^2/4;

dfdx=(zh(x(1)+eps,x(2),x(3))-zh(x(1)-eps,x(2),x(3)))/eps/2;
dfdy=(zh(x(1),x(2)+eps,x(3))-zh(x(1),x(2)-eps,x(3)))/eps/2;
dfdb=(zh(x(1),x(2),x(3)+eps)-zh(x(1),x(2),x(3)-eps))/eps/2;

jac=[d2fdx2,d2fdxdy,d2fdxdb;
    d2fdxdy, d2fdy2,d2fdydb;
    dfdx,dfdy,dfdb];

% g=[(zh(x(1)+eps,x(2))-zh(x(1)-eps,x(2)))/2/eps;(zh(x(1),x(2)+eps)-zh(x(1),x(2)-eps))/2/eps];

g=[dfdx;dfdy;zh(x(1),x(2),x(3))];

x=x-jac\g;
xa=[xa,x];
det(jac);
norm(g);
fprintf("ar: %4.2e\tai: %4.2e\tR: %4.2f\tnorm(g): %4.2e\tomi: %4.2e\n",x(1),x(2),x(3),norm(g),g(end));
%%
om=imagOmega(x(1),x(2),bbar,x(3),data);
reh=a.x(end).^2;
save("abs-reh"+num2str(reh)+".mat",'data','x','om','reh');
%% 3dof beta

clc;
% x=[0.25;-0.05];
% x=[0.2;0.2];
x=[0.2;-0.1;0.1];
xa=[x];
eps=1e-6;

zh=@(x,y,bbar) imagOmega(x,y,bbar,R);
%%


d2fdx2=(zh(x(1)+eps,x(2),x(3))-2*zh(x(1),x(2),x(3))+zh(x(1)-eps,x(2),x(3)))/eps^2;
d2fdxdy=(zh(x(1)+eps,x(2)+eps,x(3))+zh(x(1)-eps,x(2)-eps,x(3))-zh(x(1)+eps,x(2)-eps,x(3))-zh(x(1)-eps,x(2)+eps,x(3)))/eps^2/4;
d2fdy2=(zh(x(1),x(2)+eps,x(3))-2*zh(x(1),x(2),x(3))+zh(x(1),x(2)-eps,x(3)))/eps^2;

d2fdxdb=(zh(x(1)+eps,x(2),x(3)+eps)+zh(x(1)-eps,x(2),x(3)-eps)-zh(x(1)+eps,x(2),x(3)-eps)-zh(x(1)-eps,x(2),x(3)+eps))/eps^2/4;
d2fdydb=(zh(x(1),x(2)+eps,x(3)+eps)+zh(x(1),x(2)-eps,x(3)-eps)-zh(x(1),x(2)+eps,x(3)-eps)-zh(x(1),x(2)-eps,x(3)+eps))/eps^2/4;

dfdx=(zh(x(1)+eps,x(2),x(3))-zh(x(1)-eps,x(2),x(3)))/eps/2;
dfdy=(zh(x(1),x(2)+eps,x(3))-zh(x(1),x(2)-eps,x(3)))/eps/2;
dfdb=(zh(x(1),x(2),x(3)+eps)-zh(x(1),x(2),x(3)-eps))/eps/2;

jac=[d2fdx2,d2fdxdy,d2fdxdb;
    d2fdxdy, d2fdy2,d2fdydb;
    dfdx,dfdy,dfdb];

% g=[(zh(x(1)+eps,x(2))-zh(x(1)-eps,x(2)))/2/eps;(zh(x(1),x(2)+eps)-zh(x(1),x(2)-eps))/2/eps];

g=[dfdx;dfdy;zh(x(1),x(2),x(3))];

x=x-jac\g;
xa=[xa,x];
det(jac)
norm(g)


%%
hold on;
plot(xa(1,:),xa(2,:),'-x')
title("imag ev last point "+num2str(zh(xa(1,end),xa(2,end)),'%4.2e'));

%% crit points of Ray
plot(x,real(U(2:4:end)*ev(iev)+bbar*U(3:4:end)-omega)); hold on;  grid on;
plot(x,imag(U(2:4:end)*ev(iev)+bbar*U(3:4:end)-omega)); hold on;  grid on;
%%
eva=[];
% oma=oma2; lam1=0.25-0.5i; evaa=[lam1];
for i=1:length(oma)
%     for i=40:length(oma)
omega=oma(i);
[g,jac0,jac1,jac2]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
% tic; [evc,evs]=polyeig(jac0,jac1,jac2); ev=(evs); toc; 
% tic; 
[evc,evs]=polyeigs2(jac0,jac1,jac2,4,a0r+a0i*1i); ev=diag((evs));
% toc; 
eva=[eva,ev];
i
end
%%
clf;
plot(eva,'.'); %hold on;
grid on;
