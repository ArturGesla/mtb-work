clc; clear;

np=100;
x=linspace(0,20,np);
u=zeros(np*3,1);

a=load("../vonkarman-bvp4c/vk.mat"); 
u(1:3:end)=interp1(a.sol.x,a.sol.y(1,:),x);
u(2:3:end)=interp1(a.sol.x,a.sol.y(2,:),x);
u(3:3:end)=interp1(a.sol.x,a.sol.y(3,:),x);

[g,jac]=evalJacRhs(u,x);

