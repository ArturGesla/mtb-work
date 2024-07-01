   %vonK code 
   %uvw at cell wall, p at cell center, mom eqs written at wall,
   %conti at center

   % z^    ----- u,v,w | x,y mom
   %  |    | . |       | continuity     
   %  |    -----
   %  |    | . |
   %  |    -----
   %  |    | . |
   %  |-----------------> 




function [g,jac]=evalJacRhs(u,z)

np=length(z);
ii=[]; jj=[]; vv=[];
g=u*0;

%conti
for i=2:np
    iif=(i-1)*3+1; iig=(i-1)*3+2; iih=(i-1)*3+3;
    iifzm=iif-3; iihzm=(i-1)*3+3-3;
    dz=z(i)-z(i-1);
    g(iif)=2*(u(iif)+u(iifzm))/2+(u(iih)-u(iihzm))/dz;

   ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1];
   ii=[ii;iif]; jj=[jj;iifzm]; vv=[vv;1];
   ii=[ii;iif]; jj=[jj;iih]; vv=[vv;1/dz];
   ii=[ii;iif]; jj=[jj;iihzm]; vv=[vv;-1/dz];
end

%xmom ymom
for i=2:np-1
    iif=(i-1)*3+1; iig=(i-1)*3+2; iih=(i-1)*3+3;
    iifzm=iif-3; iigzm=iig-3; iihzm=iih-3; 
    iifzp=iif+3; iigzp=iig+3; iihzp=iih+3; 
    
    zp=z(i); zn=(z(i)+z(i+1))/2; zs=(z(i)+z(i-1))/2;
    

    %xmom
    fn=(u(iif)*(z(i+1)-zn)+u(iifzp)*(zn-z(i)))/(z(i+1)-z(i));
    fs=(u(iifzm)*(z(i)-zs)+u(iif)*(zs-z(i-1)))/(z(i)-z(i-1));

    dfn=(u(iifzp)-u(iif))/(z(i+1)-z(i));
    dfs=(u(iif)-u(iifzm))/(z(i)-z(i-1));

    g(iig)=u(iif)^2-u(iig)^2+(fn-fs)/(zn-zs)*u(iih)-(dfn-dfs)/(zn-zs);
   
    % 
    ii=[ii;iig]; jj=[jj;iif]; vv=[vv;2*u(iif)];
    ii=[ii;iig]; jj=[jj;iif]; vv=[vv;((z(i+1)-zn))/(z(i+1)-z(i))/(zn-zs)*u(iih)];
    ii=[ii;iig]; jj=[jj;iif]; vv=[vv;-((zs-z(i-1)))/(z(i)-z(i-1))/(zn-zs)*u(iih)];
    ii=[ii;iig]; jj=[jj;iif]; vv=[vv;(-1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iif]; vv=[vv;(1)/(z(i)-z(i-1))*(1)/(zn-zs)];

    ii=[ii;iig]; jj=[jj;iifzp]; vv=[vv;((zn-z(i)))/(z(i+1)-z(i))/(zn-zs)*u(iih)];
    ii=[ii;iig]; jj=[jj;iifzp]; vv=[vv;(1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
    
    ii=[ii;iig]; jj=[jj;iifzm]; vv=[vv;-((z(i)-zs))/(z(i)-z(i-1))/(zn-zs)*u(iih)];
    ii=[ii;iig]; jj=[jj;iifzm]; vv=[vv;(1)/(z(i)-z(i-1))*-(1)/(zn-zs)];

    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;-2*u(iig)];
    ii=[ii;iig]; jj=[jj;iih]; vv=[vv;(fn-fs)/(zn-zs)];

    %ymom
    gn=(u(iig)*(z(i+1)-zn)+u(iigzp)*(zn-z(i)))/(z(i+1)-z(i));
    gs=(u(iigzm)*(z(i)-zs)+u(iig)*(zs-z(i-1)))/(z(i)-z(i-1));

    dgn=(u(iigzp)-u(iig))/(z(i+1)-z(i));
    dgs=(u(iig)-u(iigzm))/(z(i)-z(i-1));

    g(iih)=2*u(iif)*u(iig)+(gn-gs)/(zn-zs)*u(iih)-(dgn-dgs)/(zn-zs);
    
    ii=[ii;iih]; jj=[jj;iig]; vv=[vv;2*u(iif)];
    ii=[ii;iih]; jj=[jj;iig]; vv=[vv;((z(i+1)-zn))/(z(i+1)-z(i))/(zn-zs)*u(iih)];
    ii=[ii;iih]; jj=[jj;iig]; vv=[vv;-((zs-z(i-1)))/(z(i)-z(i-1))/(zn-zs)*u(iih)];
    ii=[ii;iih]; jj=[jj;iig]; vv=[vv;(-1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iig]; vv=[vv;(1)/(z(i)-z(i-1))*(1)/(zn-zs)];

    ii=[ii;iih]; jj=[jj;iigzp]; vv=[vv;((zn-z(i)))/(z(i+1)-z(i))/(zn-zs)*u(iih)];
    ii=[ii;iih]; jj=[jj;iigzp]; vv=[vv;(1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
    
    ii=[ii;iih]; jj=[jj;iigzm]; vv=[vv;-((z(i)-zs))/(z(i)-z(i-1))/(zn-zs)*u(iih)];
    ii=[ii;iih]; jj=[jj;iigzm]; vv=[vv;(1)/(z(i)-z(i-1))*-(1)/(zn-zs)];

    ii=[ii;iih]; jj=[jj;iif]; vv=[vv;2*u(iig)];
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;(gn-gs)/(zn-zs)];



  
end

%bc
i=1; iif=(i-1)*3+1; iig=(i-1)*3+2; iih=(i-1)*3+3;
g(iif)=u(iif); g(iig)=u(iig)-1; g(iih)=u(iih);

ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1];
ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1];
ii=[ii;iih]; jj=[jj;iih]; vv=[vv;1];
  

i=np; iif=(i-1)*3+1; iig=(i-1)*3+2; iih=(i-1)*3+3;
g(iig)=u(iif); g(iih)=u(iig);

ii=[ii;iig]; jj=[jj;iif]; vv=[vv;1];
ii=[ii;iih]; jj=[jj;iig]; vv=[vv;1];


jac=sparse(ii,jj,vv);
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