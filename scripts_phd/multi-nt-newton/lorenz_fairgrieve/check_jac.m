uorg=u;

eps=1e-6;

% i=2;
up=u; um=u;
up(i)=u(i)+eps;
um(i)=u(i)-eps;

u=up;
evalJacRhs;
 gp=g;

 u=um;
evalJacRhs;
 gm=g;

evalJacRhs;
u=uorg;

col=(gp-gm)/2/eps

J(:,i)
norm(col-J(:,i))
i=i+1