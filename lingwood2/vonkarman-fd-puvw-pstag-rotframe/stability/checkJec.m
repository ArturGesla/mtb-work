np=30;
x=linspace(0,20,np);
u=rand(np*4,1);
eps=1e-3;

a=load("../vk-np-201.mat");
x=a.x;
u=a.u*0;
U=a.u;
omega=0+0.04i;
% beta=bbar=beta/R;
bbar=0.126;
R=1;
alpha=2; 



clc;
for i=1:length(u)

    up=u; up(i)=up(i)+eps;
    um=u; um(i)=um(i)-eps;

    [g,jac0,jac1]=evalJacRhsStab(u,x,U,omega,bbar,R,alpha);
    jac=jac0+alpha*jac1;
    [gp,jac0,jac1]=evalJacRhsStab(up,x,U,omega,bbar,R,alpha);
    [gm,jac0,jac1]=evalJacRhsStab(um,x,U,omega,bbar,R,alpha);
    
    
%     [gp,jacp]=evalJacRhs(up,x);
%     [gm,jacm]=evalJacRhs(um,x);

    err=(gp-gm)/2/eps-full(jac(:,i));
    fprintf("%d i \t err: %4.2e \n",i,norm(err));
    if(norm(err)>1e-10) break; end;


end