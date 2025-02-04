a=[-10 5 2 6 7]; %coeffs a0 .. an
a=[9600 -12400 2904 -105 1]; %coeffs a0 .. an x^4 - 105 x^3 + 2904 x^2 - 12400 x + 9600
an=a(end); 
n=length(a);

A=diag(ones(n-2,1),+1);
A(end,:)=-a(1:end-1)/an;
eig(A)
roots(fliplr(a))