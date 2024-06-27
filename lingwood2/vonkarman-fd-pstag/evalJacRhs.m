function [g,jac]=evalJacRhs(u,z)

np=length(z);


%conti
for i=2:np
    iif=(i-1)*3+1; iig=(i-1)*3+2; iih=(i-1)*3+3;
    iifzm=iif-3; iihzm=(i-1)*3+3-3;
    dz=z(i)-z(i-1);
    g(iif)=2*(u(iif)+u(iifzm))/2+(u(iih)-u(iihzm))/dz;
end

%xmom ymom
for i=2:np-1
    iif=(i-1)*3+1; iig=(i-1)*3+2; iih=(i-1)*3+3;
    iifzm=iif-3; iihzm=(i-1)*3+3-3;
    iifzp=iif+3; 
    
    dz=z(i)-z(i-1);
    g(iig)=u(iif)^2-u(iig)^2+(u(iifzp)-u(iifzm))/2/dz;
end

jac=0;
end

% function dydx = bvpfcn(x,y) % equation to solve
% dydx = zeros(5,1);
% dydx = [y(4)
%        y(5)
%        -2*y(1)
%        y(1)^2-y(2)^2+y(4)*y(3)
%        2*y(1)*y(2)+y(5)*y(3)];
% end
% %--------------------------------
% function res = bcfcn(ya,yb) % boundary conditions
% res = [ya(1)
%     ya(2)-1
%     ya(3)
%     yb(1)
%     yb(2) ];
% end