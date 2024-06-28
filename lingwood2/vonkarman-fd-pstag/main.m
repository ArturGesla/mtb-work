clc; clear;

np=100;
x=linspace(0,20,np);
u=zeros(np*3,1);

a=load("../vonkarman-bvp4c/vk.mat");
u(1:3:end)=interp1(a.sol.x,a.sol.y(1,:),x);
u(2:3:end)=interp1(a.sol.x,a.sol.y(2,:),x);
u(3:3:end)=interp1(a.sol.x,a.sol.y(3,:),x);

[g,jac]=evalJacRhs(u,x);

norm(g)
%%
u=rand(np*3,1);

%%
for i=1:15
    [g,jac]=evalJacRhs(u,x);
    fprintf("%d i \t norm(rhs): %4.2e \n",i,norm(g));
    u=u-jac\g;
    if(norm(g)<1e-14) break; end;

end

%%
plot(x,reshape(u,[3,np])')