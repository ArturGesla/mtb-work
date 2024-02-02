i=1;
%%

uorg=u;



eps=1e-6;

% i=2;
up=u; um=u;
up(i)=u(i)+eps;
um(i)=u(i)-eps;

u=up;
[g,jac]=calculateRhsAndJac(3,2,u);
 gp=g;

 u=um;
[g,jac]=calculateRhsAndJac(3,2,u);
 gm=g;

[g,jac]=calculateRhsAndJac(3,2,u);
J=full(jac);
u=uorg;

col=(gp-gm)/2/eps;
col=col';

J(:,i);
norm(col-J(:,i))
i=i+1