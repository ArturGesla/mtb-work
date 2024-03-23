% Renumbering
% 
% for ieq = 1:nt*nvar*nx*ny
%     [ix,jy,nv,k] = ieq2ind(ieq,nvar,nx,nt);
%     disp( [ieq,nv,k,ix,jy] )
% end

function [ix,jy,nv,k] = ieq2ind(ieq,nvar,nx,nt)
%IEQ2IND Convert ieq to indices ix,jy, k and nv variable number 
%   
nv  = mod(ieq-1,nvar)+1;
k   = mod((ieq - 1 - (nv-1))/nvar,nt)+1;
ix  = mod((ieq - 1 - (nv -1) - (k-1)*nvar)/(nvar*nt),nx)+1;
jy  = (ieq - nv - (k-1)*nvar - (ix-1)*nvar*nt)/(nt*nx*nvar) + 1;
end
