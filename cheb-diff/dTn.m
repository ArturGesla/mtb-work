function f=dTn(x,n)
% if (abs(x)~=1)
    f=n*Un(x,n-1);
% elseif(x==1)
% %     f=(n^4-n^2)/3;
% elseif(x==-1)
% %     f=(-1)^n*(n^4-n^2)/3;
% 
% end

%     f=n*((n+1)*Tn(x,n)-Un(x,n))/(x^2-1);
%     f=n*((n+1)*cos(n*acos(x))-sin((n+1)*acos(x))./sin(acos(x)))./(x.^2-1);


end