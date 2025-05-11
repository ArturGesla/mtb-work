function [evc,evs]=polyeigs2(jac0,jac1,jac2,k,sigma)

n=length(jac0);
A=[zeros(n,n), eye(n);
    jac0,jac1];
B=[-eye(n),zeros(n,n);
    zeros(n,n),jac2];

[c1,c2]=eigs((A),-(B),k,sigma);
evc=c1(1:n,:);
evs=c2;

end