% np=30;
% x=linspace(0,20,np);
% u=rand(np*4,1);
eps=1e-3;

a=load("../vk-np-50.mat");
% a=load("../vk-np-140.mat");
x=a.x;
u=a.u*0;
U=a.u;
omega=0.004i;
% beta=bbar=beta/R;
R=515; bbar=0.0117;

alpha=2+0.1i; 



clc;
for i=1:length(u)

    up=u; up(i)=up(i)+eps;
    um=u; um(i)=um(i)-eps;

    [g,jac0,jac1,jac2]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
    jac=jac0+alpha*jac1+alpha.^2*jac2;
    [gp,~,~,~]=evalJacRhsStab(up,x,U,omega,bbar,R,alpha);
    [gm,~,~,~]=evalJacRhsStab(um,x,U,omega,bbar,R,alpha);
    
    
%     [gp,jacp]=evalJacRhs(up,x);
%     [gm,jacm]=evalJacRhs(um,x);

    err=(gp-gm)/2/eps-full(jac(:,i));
    fprintf("%d i \t err: %4.2e \n",i,norm(err));
    if(norm(err)>1e-10), break; end


end
rank(full(jac))