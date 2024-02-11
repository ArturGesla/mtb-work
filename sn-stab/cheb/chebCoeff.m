function a=chebCoeff(x,f,n)
if(x(1)~=-1 ) error("err"); end
if(x(end)~=1 ) error("err"); end
% n=1;
% dx=0.001; x=0+dx:dx:pi-dx;
% x=linspace(0,pi,length(f));
x=x(2:end-1); f=f(2:end-1);
% x=fliplr(cos(x));
th=acos(x); Tn=cos(n*th);
a=trapz(x,f.*Tn./sqrt(1-x.^2))/pi*2;
if(n==0) a=a/2; end
end