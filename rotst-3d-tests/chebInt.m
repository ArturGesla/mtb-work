function I=chebInt(n)
if mod(n,2)==1
    I=0;
else
    I= 2/(1-n^2);
end