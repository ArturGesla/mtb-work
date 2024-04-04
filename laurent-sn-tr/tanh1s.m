function [dx] = tanh1s(L,n,dxl)
%TANH1S tanh distribution of step size with one prescribed step.
%     DX = TANH1S(L,N,DXL) returns a stepsize vector DX
%     for a segment of length L with N+1 grid points.
%     DXL is the first (left) prescribed step.

% 26/10/2019.  witko@limsi.fr 
% Strongly inspired from meshgen developped by Yann Fraigneau

xreg = linspace(0,1,n+1);

sn1 = dxl/L;
b = L/n/dxl;

% Try to find a good guess for alpha before solving it.

if b > 2.7829681
    v     = log(b);
    w     = 1./b - 0.02852731;
    alphaguess   = v+(1.+1./v)*log(2.*v)-0.02041793+0.24902722*w ...
        +1.9496443*w^2-2.6294547*w^3+8.56795911*w^4 ;
elseif (b >= 1.0001) && (b <= 2.7829681)
    bbar = b-1;
    alphaguess =  sqrt(6.*(bbar))*(1.-0.15*(bbar)+0.057321429*(bbar)^2 ...
        -0.024907295*(bbar)^3+0.0077424461*(bbar)^4 ...
        -0.0010794123*(bbar)^5 );
elseif (b >=0.9999 ) && (b < 1.0001)
    alphaguess = b - 1;
elseif (b >=0.26938972 ) && (b < 0.9999)
    bbar  = 1.-b;
    alphaguess   = sqrt(6.*(bbar))*(1.+0.15*(bbar)+0.057321429*(bbar)^2 ...
        +0.048774238*(bbar)^3-0.0053337753*(bbar)^4 ...
        +0.075845134*(bbar)^5 );
elseif (b < 0.26938972 )
    alphaguess   = pi*(1.-b+b^2-(1.+pi^2/6.)*b^3+6.794732*b^4 ...
        -13.205501*b^5+11.726095*b^6);
else
    error('tanh1s :  problem with b')
end

% Writing the appropriate function according to b
r = 1/n;
if (b>=1.0001)
    myfun = @(alpha) (1 + tanh(0.5*alpha*(r-1))/tanh(0.5*alpha)) - sn1;
elseif (b < 0.9999)
    myfun = @(alpha) (1 + tan(0.5*alpha*(r-1))/tan(0.5*alpha)) - sn1;
else
    myfun = @(alpha) (r*(1.0-0.5*alpha*(1.0-r)*(2.0-r))) - sn1;
end

% solve to find alpha
alpha = fzero(myfun,alphaguess);

%%
if (b>=1.0001)
    xfinal = 1 + tanh(0.5*alpha*(xreg-1))/tanh(0.5*alpha);
elseif (b < 0.9999)
    xfinal = 1 + tan(0.5*alpha*(xreg-1))/tan(0.5*alpha);
else
    xfinal =  (xreg.*(1.0-0.5*alpha*(1.0-xreg).*(2.0-xreg)));
end

xfinal = L*xfinal;

dx = diff(xfinal);

end
