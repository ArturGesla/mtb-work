function c=ci(n,m) 
if (n==m)
c=0;
else
c=(1-(-1)^(abs(n-m)))/pi/(n-m);
end
end
