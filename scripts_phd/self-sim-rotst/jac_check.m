%jacobian check
u=rand(nz*4+1,1);
[rhs,jac]=calculateJacAndRhs(zc,zw,u,Re);
eps=1e-6;
jacNum=[];

for i=1:length(u)
    u1=u; u1(i)=u1(i)+eps;
    u2=u; u2(i)=u2(i)-eps;
    
    [rhs1,jac1]=calculateJacAndRhs(zc,zw,u1,Re);
    [rhs2,jac2]=calculateJacAndRhs(zc,zw,u2,Re);
    
    col=(rhs1-rhs2)/2/eps;
    jacNum(:,i)=col;
end

%%
mesh(full(jac)-jacNum)
%%
(full(jac)-jacNum)