

function Eqk = sspconv(u,v)
% We assume u and v have the same size and are vectors 
% Here, we return only the 
% central part of the convolution of the same size than u 

% Ugly coding with loops  but OK because it is used only in
% preprocessing. 

nt = 2*length(u)-1; %  nt is here different than in the rest of the code!
nshift = (nt+1)/2; 
syms Eqk [1 length(u)];
ktab = -(nt-1)/2:(nt-1)/2;

uf = [conj(fliplr(u(2:end))) u ]; 
vf = [conj(fliplr(v(2:end))) v ];

for knum= 0:(nt-1)/2 
    kind = knum + 1;
    Eqk(kind) = 0;
    for pnum = ktab
        for lnum = ktab
            if knum==pnum+lnum
                Eqk(kind) = Eqk(kind) + uf(pnum+nshift)*vf(lnum+nshift);
            end
        end
    end
end

return
end
