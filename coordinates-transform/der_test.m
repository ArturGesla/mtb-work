clc; clear;
f=@(r) r^2;

eps=1e-3;
df=@(r) (f(r+eps)-f(r-eps))/2/eps;
d2f=@(r) (f(r+eps)+f(r-eps)-2*f(r))/eps/eps;

z=1; L=2;
fth=@(th) (L-z)*tan(th);
dfth=@(th) (fth(th+eps)-fth(th-eps))/2/eps*cos(th)^2/(L-z);
d2fth=@(th) -4*sin(th)*cos(th)/(L-z)*(fth(th+eps)-fth(th-eps))/2/eps*cos(th)^2/(L-z)+(cos(th)^2/(L-z))^2*(fth(th+eps)+fth(th-eps)-2*fth(th))/eps/eps;


% d2fth=

tha=[0:0.01:5];
va=tha*0;
for i=1:length(tha)
    va(i)=d2fth(tha(i));
end
%
plot(tha,va)