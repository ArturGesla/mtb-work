function x = besselzeroj(n,q,opt)

% first q roots of bessel function Jn(x), integer order.
% if opt = 'd', first q roots of dJn(x)/dx, integer order. 
% if opt is not provided, the default is Jn(x).
% includes the root x = 0 for dJ0(x)/dx (zeroth order)
% but not for higher orders. 
%
% starting point for for zeros of Jn is borrowed from Cleve Moler, 
% but starting points for both Jn and Jn' can be found in
% Abramowitz and Stegun 9.5.12, 9.5.13. 
%
% David Goodmanson
%
% function x = besselzeroj(n,q,opt)

k = 1:q;
if nargin==3 & opt=='d'
  % find starting points for J'
  beta = (k + n/2 - 3/4)*pi;
  mu = 4*n^2;
  x = beta - (mu+3)./(8*beta) - 4*(7*mu^2+82*mu-9)./(3*(8*beta).^3);
  % Newton
  for j=1:8
    xnew = x - besseljdv(n,x)./ ...
        (besselj(n,x).*((n^2./x.^2)-1) -besseljdv(n,x)./x);
    x = xnew;    
  end
else
  % find starting points for J
  beta = (k + n/2 - 1/4)*pi;
  mu = 4*n^2;
  x = beta - (mu-1)./(8*beta) - 4*(mu-1)*(7*mu-31)./(3*(8*beta).^3);
  % Newton
  for j=1:8
    xnew = x - besselj(n,x)./besseljdv(n,x);
    x = xnew;
  end
end
end

function y = besseljdv(n,x)

% derivative of bessel function of integer order,
% plain vanilla version 
% 
% function y = besseljdv(n,x)

y = -besselj(n+1,x) + n*besselj(n,x)./x;

% get rid of nans, integer case
if n==1
  y(x==0) = 1/2;
else
  y(x==0) = 0;
end
end
