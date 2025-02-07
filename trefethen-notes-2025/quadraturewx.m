function [w,x]=quadraturewx(N)
beta=0.5./sqrt(1-(2*(1:N)).^(-2));
T=diag(beta,1)+diag(beta,-1); % rec relation coeffs
[V,D]=eig(T);
x=diag(D); [x,i]=sort(x);
w=2*V(1,i).^2;
end