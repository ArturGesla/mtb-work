i=1;
%%

uorg=u; 
mu=0.04; gmm=1;
u=rand(nt*3*2+1,1);


eps=1e-6;

% i=2;
up=u; um=u;
up(i)=u(i)+eps;
um(i)=u(i)-eps;

u=up;
[g,jac]=calculateRhsAndJac(3,nt,u,mu,gmm);
 gp=g;

 u=um;
[g,jac]=calculateRhsAndJac(3,nt,u,mu,gmm);
 gm=g;

[g,jac]=calculateRhsAndJac(3,nt,u,mu,gmm);
J=full(jac);
u=uorg;

col=(gp-gm)/2/eps;
col=col';

J(:,i);
i
norm(col-J(:,i))
i=i+1;