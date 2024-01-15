i=1;
%%

uorg=u;



eps=1e-6;

% i=2;
up=u; um=u;
up(i)=u(i)+eps;
um(i)=u(i)-eps;

u=up;
[g,jac]=calculateRhsAndJac(3,nt,u);
 gp=g;

 u=um;
[g,jac]=calculateRhsAndJac(3,nt,u);
 gm=g;

[g,jac]=calculateRhsAndJac(3,nt,u);
J=full(jac);
u=uorg;

col=(gp-gm)/2/eps;
col=col';

J(:,i);
i
norm(col-J(:,i))
i=i+1;