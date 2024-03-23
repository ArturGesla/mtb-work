function Eqk = sspconv(u,v)
% We assume u and v have the same size and are vectors 
% Here, we return only the 
% central part of the convolution of the same size than u 

% Ugly coding with loops  but OK for now

nt = length(u);
if mod(nt,2) == 0
    error(' Error nt must be odd')
end
nshift = (nt+1)/2; 
syms Eqk [1 nt];
ktab = -(nt-1)/2:(nt-1)/2;

for knum= ktab
    kind = knum + nshift;
    Eqk(kind) = 0;
    for pnum = ktab
        for lnum = ktab
            if knum==pnum+lnum
                Eqk(kind) = Eqk(kind) + u(pnum+nshift)*v(lnum+nshift);
            end
        end
    end
end
% disp('Warning ssp')
% Eqk = u.*v;
% 
return
end