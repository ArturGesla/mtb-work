function [ix,jy,nv,k] = ieq2ind(ieq,nvar,nx,nt)
%IEQ2IND Convert ieq to indices ix,jy, k and nv variable number
%
nv  = mod(ieq-1,nvar)+1;
k   = mod((ieq - 1 - (nv-1))/nvar,2*nt)+1;
ix  = mod((ieq - 1 - (nv -1) - (k-1)*nvar)/(nvar*2*nt),nx)+1;
jy  = (ieq - nv - (k-1)*nvar - (ix-1)*nvar*2*nt)/(2*nt*nx*nvar) + 1;
end

