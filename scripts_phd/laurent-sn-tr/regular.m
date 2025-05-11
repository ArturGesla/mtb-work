function [dx] = regular(L,n)
%REGULAR returns a regular distribution of step size.
%     DX = REGULAR(L,N) returns a N vector DX with step size
%     L/N

% 26/10/2019.  witko@limsi.fr 


xreg = linspace(0,L,n+1);
dx = diff(xreg);


end