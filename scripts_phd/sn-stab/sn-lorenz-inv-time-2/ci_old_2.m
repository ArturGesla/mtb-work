function c=ci(n,m) 
nt=7;
if (n==m)
c=0;
elseif(abs(n)>(nt-1))
c=0;
else
c=(1-(-1)^(abs(n-m)))/pi/(n-m);    
end

end
