function [dx] = mixedCosReg(L,n,beta)
%REGULAR returns a regular distribution of step size.
%     DX = REGULAR(L,N) returns a N vector DX with step size
%     L/N

% 26/10/2019.  witko@limsi.fr 


xreg = linspace(0,L,n+1);
xcos=(-cos([linspace(0,pi,n+1)])+1)/2*L;
xreg;
xcos;
dx = diff(xreg*(1-beta)+xcos*(beta));


end
