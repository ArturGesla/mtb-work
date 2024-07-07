function [evc,evs]=polyeigs(jac0,jac1,jac2,jac3,k,sigma)

n=length(jac0);
A=[zeros(n,n), eye(n),zeros(n,n);
    zeros(n,n),zeros(n,n),eye(n);
    jac0,jac1,jac2];
B=[-eye(n),zeros(n,n),zeros(n,n);
    zeros(n,n),-eye(n),zeros(n,n);
    zeros(n,n),zeros(n,n),jac3];

[c1,c2]=eigs((A),-(B),k,sigma);
evc=c1(1:n,:);
evs=c2;

end