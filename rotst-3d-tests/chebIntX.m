function I=chebIntX(n)
if mod(n,2)==0
    I=0;
elseif n==1 
    I= 2/3;
else
    I= 2/(1-n^2)-(1/(n+1)*chebInt(n+1)-1/(n-1)*chebInt(n-1))/2;
end