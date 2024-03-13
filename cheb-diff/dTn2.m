function f=dTn2(x,n)
if (length(x)==1)
f=n*((n+1)*Tn(x,n)-Un(x,n))/(x^2-1);
end