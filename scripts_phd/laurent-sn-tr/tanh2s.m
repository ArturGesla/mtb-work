function [dx] = tanh2s(L,n,dxl,dxr)
%TANH2S tanh distribution of step size with two prescribed steps.
%     DX = TANH1S(L,N,DXL) returns a stepsize vector DX
%     for a segment of length L with N+1 grid points
%     DXL is the first (left) prescribed step.
%     DXR is the last (right) prescribed step.

% 26/10/2019.  witko@limsi.fr 
% Strongly inspired from meshgen developped by Yann Fraigneau

%L = 0.071428571428571 ;
%n = 128;
%dxl = 3.5e-4;
%dxr = 4.e-4;


xreg = linspace(0,1,n+1);

s1    = L/dxl/n;
s2    = L/dxr/n;
s1n   = dxl/L;
s2n   = dxr/L;

a = sqrt(s1/s2);
b = sqrt(s1*s2);


% Try to find good guesses before solving with Newton.


if(b >= 1.0001)  %-- first case > 1.0001
    
    if(b < 2.7829681)
        
        %-- initial guess by series expansion
        
        bbar  = b-1;
        delta = sqrt(6.*(bbar))*(1.-0.15*(bbar)+0.057321429*(bbar)^2 ...
        -0.024907295*(bbar)^3+0.0077424461*(bbar)^4 ...
        -0.0010794123*(bbar)^5 );
    else
        
        %-- initial guess by series expansion
        
        v = log(b);
        w = 1./b - 0.02852731;
        delta = v+(1.+1./v)*log(2.*v)-0.02041793+0.24902722*w ...
        +1.9496443*w^2-2.6294547*w^3+8.56795911*w^4;
        
        
    end
    
elseif (b < 0.9999) %-- second  case b < 0.9999
    
    if(b < 0.26938972)
        
        %-- initial guess by series expansion
        
        pi    = acos(0.)*2.;
        delta = pi*(1.-b+b^2-(1.+pi^2/6.)*b^3+6.794732*b^4 ...
        -13.205501*b^5+11.726095*b^6);
    else
        
        %-- initial guess by series expansion
        
        bbar  = 1.-b;
        delta = sqrt(6.*(bbar))*(1.+0.15*(bbar)+0.057321429*(bbar)^2 ...
        +0.048774238*(bbar)^3-0.0053337753*(bbar)^4 ...
        +0.075845134*(bbar)^5 );
    end
    
else %--- third case 0.9999<=b<=1.0001
    
    %-- initial guess by series expansion
    
    delta = b-1.0;
end

%-- find root delta and a to obtain exactly s1 and s2
%-- use of Newton-Raphson method for zero search

x1    = a;
x2    = delta;
f1    = 0.0;
f2    = 0.0;
a11   = 1.0;
a12   = 0.0;
a21   = 0.0;
a22   = 1.0;

itmax  = 80; % there is something strange that so many iterations are needed
itmax1 = 20;

for it = 1:itmax
    
    % Compute Functions and Jacobian
    r1  =  1/n;
    r2  = (n-1)/n;
    
    if (b >= 1.0001)
        
        u1 = 0.5*(1.0+tanh(x2*(r1-0.5))/tanh(0.5*x2));
        u2 = 0.5*(1.0+tanh(x2*(r2-0.5))/tanh(0.5*x2));
        f1 = (x1 + (1.0-x1)*u1)*s1n - u1;
        f2 = (x1 + (1.0-x1)*u2)*s2n - x1*(1.0-u2);
        du1dx2 = 0.5*(r1-0.5)/(cosh(x2*(r1-0.5)))^2/tanh(x2*0.5) ...
                - 0.25*tanh(x2*(r1-0.5))/(sinh(x2*0.5))^2;
        du2dx2 = 0.5*(r2-0.5)/(cosh(x2*(r2-0.5)))^2/tanh(x2*0.5) ...
                - 0.25*tanh(x2*(r2-0.5))/(sinh(x2*0.5))^2;
 
    elseif (b < 0.9999)
        
        u1 = 0.5*(1.0+tan(x2*(r1-0.5))/tan(0.5*x2));
        u2 = 0.5*(1.0+tan(x2*(r2-0.5))/tan(0.5*x2));
        f1 = (x1 + (1.0-x1)*u1)*s1n - u1;
        f2 = (x1 + (1.0-x1)*u2)*s2n - x1*(1.0-u2);
        du1dx2 = 0.5*(r1-0.5)/(cos(x2*(r1-0.5)))^2/tan(x2*0.5) ...
                - 0.25*tan(x2*(r1-0.5))/(sin(x2*0.5))^2;
        du2dx2 = 0.5*(r2-0.5)/(cos(x2*(r2-0.5)))^2/tan(x2*0.5) ...
                - 0.25*tan(x2*(r2-0.5))/(sin(x2*0.5))^2;

    else
        
        u1 = r1*(1.0+2.0*x2*(r1-0.5)*(1.0-r1));
        u2 = r2*(1.0+2.0*x2*(r2-0.5)*(1.0-r2));
        f1 = (x1 + (1.0-x1)*u1)*s1n - u1;
        f2 = (x1 + (1.0-x1)*u2)*s2n - x1*(1.0-u2);
        du1dx2 = 2.0*r1*(r1-0.5)*(1.0-r1);
        du2dx2 = 2.0*r2*(r2-0.5)*(1.0-r2);

    end
    a11     = (1.0-u1)*s1n;
    a12     = du1dx2*((s1n-1.0) - x1* s1n     );
    a21     = (s2n-1.0)*(1.0-u2);
    a22     = du2dx2*( s2n      - x1*(s2n-1.0));

 % solve    
    omega1= it/(itmax1-1);
    omega2= it/(itmax1-1);
    omega1= min(1.0,omega1);
    omega2= min(1.0,omega2);
    delta = (a11*a22 - a21*a12);
    rdelta= 1.0/delta;
    
    dx1   = -(f1*a22 - f2*a12)*rdelta;
    dx2   = -(f2*a11 - f1*a21)*rdelta;
    dx1   = omega1*dx1;
    dx2   = omega2*dx2;
    
    %x1old = x1;
    %x2old = x2;
    x1    = x1 + dx1;
    x2    = x2 + dx2;
    x1    = max(x1,eps);
    x2    = max(x2,eps);
    
    % if ( (abs(dx1) < eps && abs(dx2) < eps) && (abs(f1 ) < eps && abs(f2 ) < eps) ) exit
    
end

a      = x1;
delta  = x2;







% %%
if (b>=1.0001)
    u = 0.5*(1.+tanh(delta*(xreg-0.5))/tanh(0.5*delta));
elseif (b < 0.9999)
    u = 0.5*(1.+tan(delta*(xreg-0.5))/tan(0.5*delta));
else
    u = xreg.*(1.+2*delta*(xreg-0.5).*(1-xreg));
end

xfinal = L*(u./(a+(1.-a)*u));

% 
 dx = diff(xfinal);
% 
end
