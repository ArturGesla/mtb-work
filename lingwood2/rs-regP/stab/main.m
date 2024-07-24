 clc; clear;
 
 
 addpath C:\Users\Artur\Documents\GitHub\rotst2\scripts\source_for_mtb;
%%
% a=load("../vk-np-50.mat");
% a=load("../vk-np-100.mat");
% a=load("../vk-np-200.mat");
% a=load("../vk-np-400.mat");

% a=load("../vk-np-110.mat");
% a=load("../vk-np-120.mat");
a=load("../vk-np-130.mat");
% a=load("../vk-np-140.mat");
% a=load("../vk-np-180.mat");

x=a.x;
u=a.u*0;
U=a.u;
z=x;
zw=[2*z(1)-z(2),z(1:end-1),2*z(end-1)-z(end-2)];
zc=(zw(1:end-1)+zw(2:end))/2;
zw=zw(2:end);
% omega=0+0.04i;
% % beta=bbar=beta/R;
% bbar=0.126;
% R=1;
% alpha=1; 

% [g,jac0,jac1,jac2]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);

% R=515; bbar=0.0117;
% beta=0; bbar=beta/R;


%%

evaM=[];
%%
% omega=0.01i/2/2; 

% R=515; bbar=0.0117; omega=-0.025+0.004i; omega=-0.025+0.0005i; omega=-0.025-0.001i;
% R=530; beta=67; bbar=beta/R; omega=-0.025*2+0.1i/2/2/2/2/2/2/2; %omega=-0.025+0.0005i; omega=-0.025-0.001i;

R=530; beta=67; bbar=beta/R; omega=-0.025*2+0.1i/2/2/2/2/2/2/2; %omega=-0.025+0.0005i; omega=-0.025-0.001i;

ns=100; dom=0.01/2/2/4; om0r=-0.03125; om0i=4e-4; a0r=0.22; a0i=-0.14;
oma=linspace(om0r-ns/2*dom,om0r+ns/2*dom,ns)+om0i*1i;
% oma=-0.0692 :0.002/4:0.0262; oma1=oma+0.04i;
% oma=-0.0692 :0.002/4:0.0262; oma2=oma+0.015i;
% % om0=-0.0692; om1=0.0262; omx=-1:2/50:1; omx=omx.^3; omx=(omx+1)/2*(om1-om0)+om0; oma3=omx+0.0133i;
% oma=-0.0692 :0.002/4:0.0262; oma3=oma+0.0133i;
%
% ev=eig(full(jac0),-full(jac1));
 
R=27.4; bbar=0.1152; omegaa=[-0.1:0.01/4:0.1]+0.01i/2/2/4/4; alpha=0;


%%
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
evaM=[evaM,eva];

%% alplot

clf;
plot(eva,'k.'); hold on;
% plot(eva,'k.'); hold on; text(real(eva),imag(eva),num2str(im))
plot(evaM,'.'); hold on;
% plot(eva,'x'); hold on;
 plot(ev,'x'); hold on; text(real(ev),imag(ev),num2str([1:length(ev)]'))
% plot(ev,'o');
% xlim([0 0.35]);
ylim([-0.2 0.2]*2); 
grid on;
% xlim([-0.1 1.5]); ylim([-1 1]); grid on;
% xlim([-0.1 0.5]*10); ylim([-0.5 0.5]*10); grid on;
% xlim([-0.1 0.5]*40); ylim([-0.5 0.5]*40)

% pbaspect([5 1 1]); xlabel("a_r"); ylabel("a_i"); 
% % fnts=12; jfm_plt_aid_comm;
% % exportgraphics(gcf,"spatialbranches-c.eps")
% % exportgraphics(gcf,"spatialbranches-b.eps")
% exportgraphics(gcf,"spatialbranches-a.eps")

% save("eva-"+num2str(length(a.x))+".mat",'eva');

%%
clf;
iev=1;
up=reshape(real(evc(:,iev)),[4,length(x)])';
plot(up(:,[1:3]),zc,'x-'); grid on;
hold on;
plot(up(:,4),zw,'x-');
% plot(up(:,1:1),x,'x-'); grid on;
title("ev "+num2str(iev)+":"+num2str(ev(iev),'%4.2e'))
iev=iev+1
ylim([0 max(zw)])

%%
ai=0;
eva=[];
z=[];

ai=-0.5:0.1:0.2;
ar=-0.05:0.01:0.45; 


% ai=-0.2:0.01:0;
% ar=0.1:1e-2:0.3;

% R=560; bbar=0.10;
% % ai=[0];  ar=-0.05:0.01:0.45; 
% ai=[-0.2:0.01:0];  ar=0:0.01:0.4; 

R=27.4; bbar=0.1152;
ai=[0];  ar=-0.5:0.01:1.5; 
% ai=[-0.5:0.01:0.5];  ar=0:0.01:1; 



for ii=1:length(ai) 
for ir=1:length(ar) 
    alpha=ar(ir)+1i*ai(ii);
     [g,jac0,jac1,jac2,jacom]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
    jac=jac0+alpha*jac1+alpha.^2*jac2;
    [ev]=eigs(jac,-jacom,5,0.1i);
    eva=[eva,ev];
%     disp(alpha);
fprintf("ir %d\t ii %d\n",ir,ii);
    [a,b]=max(imag(ev));
    z(ir,ii)=ev(b);
    end
end

%%
save("aiarz-"+num2str(bbar)+"-"+num2str(R)+".mat",'ai','ar','z');
%%
plot(eva,'k.'); grid on;

%%
clf;
contour(ar,ai,imag(z)',40); title("Imaginary part of most unstable $\omega(\alpha)$.");
hold on; contour(ar,ai,imag(z)',[0 0],'-k');
% contour(ar,ai,real(z)',400); title("Real part of most unstable $\omega(\alpha)$.");

% surf(ar,ai,imag(z)');

xlabel("$\alpha_r$");
ylabel("$\alpha_i$");
fnts=12; jfm_plt_aid_comm; 
size_sq; colorbar(); colormap(hsv(8));

%%

exportgraphics(gcf,"saddlepoint-4.eps");


%% precise saddle

clc;
% x=[0.25;-0.05];
% x=[0.2;0.2];
x=[0.2;-0.1];
xa=[x];
eps=1e-6;
zh=@(x,y) imagOmega(x,y,bbar,R);

%%

d2fdx2=(zh(x(1)+eps,x(2))-2*zh(x(1),x(2))+zh(x(1)-eps,x(2)))/eps^2;
d2fdxdy=(zh(x(1)+eps,x(2)+eps)+zh(x(1)-eps,x(2)-eps)-zh(x(1)+eps,x(2)-eps)-zh(x(1)-eps,x(2)+eps))/eps^2/4;
d2fdy2=(zh(x(1),x(2)+eps)-2*zh(x(1),x(2))+zh(x(1),x(2)-eps))/eps^2;

jac=[d2fdx2,d2fdxdy; d2fdxdy, d2fdy2];
g=[(zh(x(1)+eps,x(2))-zh(x(1)-eps,x(2)))/2/eps;(zh(x(1),x(2)+eps)-zh(x(1),x(2)-eps))/2/eps];
x=x-jac\g;
xa=[xa,x];
det(jac)
norm(g)

%% 3dof

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


%%

clf;
a=load("boma2.mat"); plot(a.eva); hold on;
a=load("boma2-140.mat"); plot(a.eva);
a=load("boma2-180.mat"); plot(a.eva);


%%
clf; hold on;
a=load("eva-50.mat"); plot(a.eva,'g.');
a=load("eva-100.mat"); plot(a.eva,'r.');
a=load("eva-200.mat"); plot(a.eva,'b.');

%%

clf; hold on;
a=load("eva-50.mat"); plot(a.eva,'g.');
a=load("eva-100.mat"); plot(a.eva,'r.');
a=load("eva-200.mat"); plot(a.eva,'b.');
a=load("eva-110.mat"); plot(a.eva,'c.');
a=load("eva-120.mat"); plot(a.eva,'k.');