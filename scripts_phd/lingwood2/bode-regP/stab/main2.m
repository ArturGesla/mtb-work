clc; clear;
addpath C:\Users\Artur\Documents\GitHub\rotst2\scripts\source_for_mtb;
%%
a=load("../vk-np-100-k-1-L-30.mat"); a.U=a.u; z=a.x; zw=[2*z(1)-z(2),z(1:end-1),2*z(end-1)-z(end-2)]; a.zc=(zw(1:end-1)+zw(2:end))/2; a.zw=zw(2:end);
a30=a;

a=load("../vk-np-200-k-1-L-60.mat"); a.U=a.u; z=a.x; zw=[2*z(1)-z(2),z(1:end-1),2*z(end-1)-z(end-2)]; a.zc=(zw(1:end-1)+zw(2:end))/2; a.zw=zw(2:end);
a60=a;

% a=load("../vk-np-200-k-1-L-60.mat");
% a=load("../vk-np-1600-k-1-L-480.mat");

k=a.k;
up=reshape(a.u,[4,length(a.u)/4])';
plot(up(:,1:end),a.x,'k-');
%% temp branches

omega=0;
shift=0.1i+0.1;
shift=0;
ai=0;
eva=[];
z=[];

R=27.4; bbar=0.1152; beta=bbar*R;
ai=[0];  ar=-0.5:0.01:1.5; ar=1.5;
% ai=[-0.5:0.01:0.5];  ar=0:0.01:1;
R=R/sqrt(k); bbar=beta/R;ar=ar.*sqrt(k); shift=shift/k^(3/2);


for ii=1:length(ai)
    for ir=1:length(ar)
        alpha=ar(ir)+1i*ai(ii);
        
        a=a30;
        [g,jac0,jac1,jac2,jacom]=evalJacRhsStab(a.u,a.x,a.U,omega,bbar,R,alpha);
        jac=jac0+alpha*jac1+alpha.^2*jac2;
        [evc,ev]=eigs(jac,-jacom,30,shift); ev=diag(ev); %eva=[eva,ev];
        ev30=ev;
        
        a=a60;
        [g,jac0,jac1,jac2,jacom]=evalJacRhsStab(a.u,a.x,a.U,omega,bbar,R,alpha);
        jac=jac0+alpha*jac1+alpha.^2*jac2;
        [evc,ev]=eigs(jac,-jacom,30,shift); ev=diag(ev); %eva=[eva,ev];
        ev60=ev;
        
%         b=sum((abs(dev)<1e-4));
        [a2,b2]=sort(min(abs(dev)));
%         bb=b==true;
%         ev=ev(b2(1:3));
        ev=ev(b2(1:end));
        
        eva=[eva,ev];
        
        fprintf("ir %d\t ii %d\n",ir,ii);
        [a,b]=max(imag(ev));
        z(ir,ii)=ev(b);
    end
end

%%
clf;
plot(eva,'k.'); grid on; hold on;
plot(evb,'b.'); grid on; hold on;
plot(ev,'r.'); hold on; text(real(ev),imag(ev),string(num2str([1:length(ev)]')));
%%
% save('evb','evb')
%%
up=real(reshape(evc(:,30),[4,length(a60.x)]))';
plot(up,a60.x)
% semilogx(abs(up),a60.x)
