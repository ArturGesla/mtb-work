np=30;
x=linspace(0,20,np);
u=rand(np*4,1);
eps=1e-3;

clc;
tic;
for i=1:length(u)

    up=u; up(i)=up(i)+eps;
    um=u; um(i)=um(i)-eps;

    [g,jac]=evalJacRhs(u,x);
    [gp,jacp]=evalJacRhs(up,x);
    [gm,jacm]=evalJacRhs(um,x);

    err=(gp-gm)/2/eps-full(jac(:,i));
    fprintf("%d i \t err: %4.2e \n",i,norm(err));
    if(norm(err)>1e-10) break; end;


end
toc;
rank(full(jac))