
clc; clear;
%% grid
ra=1;
kr=0;
Re=1000;
L=1;
zw=(0:1/160:L)';
zw=(1-cos(zw*pi))/2;
zw=[zw; zw(end)*2-zw(end-1)]; zc=[(3*zw(1)-zw(2))/2; zw(1:end-1)+diff(zw)/2;];
nz=length(zw);
u=rand(nz*4+1,1)*0;
u(3:4:end)=-1*(u(3:4:end)+1);
% u(3:4:end)=flipud(-log(zc+0.1)).*sin(zw);
% 
%%
% Re=Re+300;
for i=1:10
tic; [rhs,jac,B]=calculateJacAndRhs(zc,zw,u,Re,ra,kr,L); %toc; 
tic; du=jac\(-rhs'); %toc;
u=u+du;
% norm(du)
% norm(rhs)
fprintf("i: %1.0d norm du: %4.2e norm rhs %4.2e\n",i,norm(du),norm(rhs))
end
%%
% uPlot=uEvReal;
% uPlot=uEvImag;
uPlot=u;
 close all;
hold on;
leg=[];
symbol="-";
plot(zc,uPlot(1:4:end-1)); leg=[leg;"p"];
plot(zc,uPlot(2:4:end-1),symbol); leg=[leg;"u"];
plot(zc,uPlot(3:4:end-1),symbol); leg=[leg;"v"];
% plot(zc,uPlot(3:4:end-1),symbol); leg=[leg;string(num2str(Re))];
plot(zw,uPlot(4:4:end-1),symbol); leg=[leg;"w"];
legend(leg,"Location","northwest"); grid on; grid minor;
title("Re: "+num2str(Re)+" | k^{1/2}: "+num2str(sqrt(u(end))));
xlim([0 1])
exportgraphics(gcf,"baseFlow-"+num2str(Re)+".png","Resolution",150);

%% evplot

% b=2;
uEvReal=real(evc(:,b));
uEvImag=imag(evc(:,b));

uPlot=uEvReal;
uPlot2=uEvImag;
% uPlot=u;
 close all;
hold on;
leg=[];
symbol="-";
plot(zc,uPlot(1:4:end-1),symbol); leg=[leg;"p"];
% plot(zc,uPlot2(1:4:end-1),symbol+"r"); leg=[leg;"p"];
plot(zc,uPlot(2:4:end),symbol); leg=[leg;"u"];
plot(zc,uPlot(3:4:end),symbol); leg=[leg;"v"];
% plot(zc,uPlot2(3:4:end),symbol+"r"); leg=[leg;"v"];
% plot(zc,uPlot(3:4:end-1),symbol); leg=[leg;string(num2str(Re))];
plot(zw,uPlot(4:4:end),symbol); leg=[leg;"w"];
legend(leg); grid on;
xlim([0 1]);

title("most unstable eigenmode | Re: "+num2str(Re)+" | ev: "+num2str(evs(b)))
exportgraphics(gcf,"eigenmode-"+num2str(Re)+".png","Resolution",150);

%%

ev=eig(full(jac),full(B));
%%
close all;
plot(ev,'x'); xlim([-20 1]); grid on;

%%
al=0.34;
rl=21.6;
%%
al=0.3;
rl=21;

k=0.313;

delta=L/sqrt(Re*k);
kr=al/delta;
ra=rl*delta;

%
% kr=6;
%  ra=1.2;
 %
[rhs,jac,B]=calculateJacAndRhs(zc,zw,u,Re,ra,kr,L); 
%
%   evs=eig(full(jac),full(B)); plot(evs,'x');
[evc,evs]=eigs((jac),(B),40,"smallestabs"); evs=diag(evs); plot(evs,'o');
% close all;

xlim([-2 1]); grid on; hold on;

% k=0.313;
% 
% delta=L/sqrt(Re*k);
% rl=ra/delta;
% display(rl)
% al=kr*delta;
% display(al)

%%
[a,b]=max(real(evs))
plot(evs(b),'+')
%%

b=2;
uEvReal=real(evc(:,b));
uEvImag=real(evc(:,b));