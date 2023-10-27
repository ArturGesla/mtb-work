i=1;
%%

uorg=u;



eps=1e-6;

% i=2;
up=u; um=u;
up(i)=u(i)+eps;
um(i)=u(i)-eps;

u=up;
% [g,jac]=calculateRhsAndJac(3,2,u);
evalJacRhs3
 gp=g;

 u=um;
% [g,jac]=calculateRhsAndJac(3,2,u);
evalJacRhs3
 gm=g;

% [g,jac]=calculateRhsAndJac(3,2,u);
evalJacRhs3
J=full(J);
u=uorg;

col=(gp-gm)/2/eps;
% col=col';
i
J(:,i);
dJ=col-J(:,i);
norm(dJ(3:3:end-1))
norm(dJ(:))
i=i+1;