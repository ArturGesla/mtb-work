function c=ci(n,m,nt) 
% function c=ci(n,m) 
% % nt=7;
% if (n==m)
% c=0;
% elseif(abs(n)>(nt-1))
% % %     elseif(abs(n)>(nt-8))
% c=0;(1-(-1)^(abs(n-m)))/pi/(n-m);    %0;
% else
% c=0;(1-(-1)^(abs(n-m)))/pi/(n-m);    
% end



if (n==m)
c=0;
elseif(abs(n)>(nt-1))
% %     elseif(abs(n)>(nt-8))
c=2*(1-(-1)^(abs(n-m)))/pi/(n-m);    %0;
else
c=2*(1-(-1)^(abs(n-m)))/pi/(n-m);    
end

% c=0;
end
