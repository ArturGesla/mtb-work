function a=get_rs(Reh,n)
% clc; clear;
%
% np=400; x=linspace(0,20,np); x=(x/20).^2*20;
%L=1; np=100+L; x=linspace(0,L,np); x=(x/L).^1*L; Re=1000; s="x"; reh=Re;
 L=(sqrt(Reh)); np=n+1; x=linspace(0,L,np); x=[x,x(end)*2-x(end-1)]; np=length(x); x=(x/L).^1*L; Re=1; s="x"; reh=L.^2;
 
% L=(sqrt(1)); np=161; x=linspace(0,L,np); x=[x,x(end)*2-x(end-1)]; np=length(x);  x=(x/L).^1*L; Re=150; s="x"; reh=Re;

% disp("Reh: "+num2str(reh));
% x=0:0.1:20; np=length(x);
% al=0.01; x=(exp(al*x)./exp(al*20)*2-1)*20;
u=zeros(np*4,1);
 u(1:4:end)=0;
u(2:4:end)=0;
u(3:4:end)=0.313;
u(end)=0.313;
% u(4:4:end)=1;


% a=load("../vonkarman-bvp4c/vk.mat");
% u(2:4:end)=interp1(a.sol.x,a.sol.y(1,:),x);
% u(3:4:end)=interp1(a.sol.x,a.sol.y(2,:),x)-1;
% u(4:4:end)=interp1(a.sol.x,a.sol.y(3,:),x);
% 
% [g,jac]=evalJacRhs(u,x);
% 
% norm(g)
%
% u=rand(np*4,1);

%
% for ii=1:1%100
for i=1:25
    [g,jac]=evalJacRhs(u,x,Re);
%     fprintf("%d i \t norm(rhs): %4.2e \t rms norm: %4.2e \n",i,norm(g),norm(g));
    u=u-jac\g;
    if(norm(g)<1e-13) break; end;

end
 if(norm(g)>1e-10) error("no conv"); end;
%  save("rs-np-"+num2str(np)+"-k-"+num2str(k)+"-reh-"+num2str(reh)+".mat",'u','np','x','k');
k=0.313;

a.u=u;
a.np=np;
a.x=x;
a.k=k;
end