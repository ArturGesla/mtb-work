   %vonK code 
   
   % z^    ----- w | z mom
   %  |    | . | p,u,v,      | continuity,xmom,ymom     
   %  |    -----
   %  |    | . |
   %  |    -----
   %  |    | . |
   %  |----------------->  
%          | . | u,v,dummy p
%          -----



function [g,jac]=evalJacRhs(u,zc,zw)

Nz=length(zc);
ii=[]; jj=[]; vv=[];
g=u*0;

%conti
for i=2:Nz-1
    iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
    iifzm=iif-4; iihzm=(i-1)*4+4-4;
    dz=zw(i)-zw(i-1);
    
    
    g(iip)=2*u(iif)+(u(iih)-u(iihzm))/dz;

   ii=[ii;iip]; jj=[jj;iif]; vv=[vv;2];
   ii=[ii;iip]; jj=[jj;iih]; vv=[vv;1/dz];
   ii=[ii;iip]; jj=[jj;iihzm]; vv=[vv;-1/dz];
end
% i=2; iip=(i-1)*4+1; g(iip)=u(iip); ii=[ii;iip]; jj=[jj;iip]; vv=[vv;1];

% %xmom 
for i=2:Nz-1
    iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
    iifzm=iif-4; iigzm=iig-4; iihzm=iih-4; 
    iifzp=iif+4; iigzp=iig+4; iihzp=iih+4; 
    iipzp=iip+4;
    
    zp=z(i); zn=(z(i)+z(i+1))/2; zs=(z(i)+z(i-1))/2;
    

    %xmom
    fn=(u(iif)*(z(i+1)-zn)+u(iifzp)*(zn-z(i)))/(z(i+1)-z(i));
    fs=(u(iifzm)*(z(i)-zs)+u(iif)*(zs-z(i-1)))/(z(i)-z(i-1));

    dfn=(u(iifzp)-u(iif))/(z(i+1)-z(i));
    dfs=(u(iif)-u(iifzm))/(z(i)-z(i-1));

    g(iif)=u(iif)^2-u(iig)^2+(fn-fs)/(zn-zs)*u(iih)-(dfn-dfs)/(zn-zs);
    % g(iig)=g(iig)-u(iig)^2-2*u(iig)-1;
    g(iif)=g(iif)-2*u(iig)-1;
   
    % 
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;2*u(iif)];
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;((z(i+1)-zn))/(z(i+1)-z(i))/(zn-zs)*u(iih)];
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;-((zs-z(i-1)))/(z(i)-z(i-1))/(zn-zs)*u(iih)];
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;(-1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;(1)/(z(i)-z(i-1))*(1)/(zn-zs)];

    ii=[ii;iif]; jj=[jj;iifzp]; vv=[vv;((zn-z(i)))/(z(i+1)-z(i))/(zn-zs)*u(iih)];
    ii=[ii;iif]; jj=[jj;iifzp]; vv=[vv;(1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
    
    ii=[ii;iif]; jj=[jj;iifzm]; vv=[vv;-((z(i)-zs))/(z(i)-z(i-1))/(zn-zs)*u(iih)];
    ii=[ii;iif]; jj=[jj;iifzm]; vv=[vv;(1)/(z(i)-z(i-1))*-(1)/(zn-zs)];

    ii=[ii;iif]; jj=[jj;iig]; vv=[vv;-2*u(iig)];
    ii=[ii;iif]; jj=[jj;iih]; vv=[vv;(fn-fs)/(zn-zs)];

    ii=[ii;iif]; jj=[jj;iig]; vv=[vv;-2];


    %ymom
    gn=(u(iig)*(z(i+1)-zn)+u(iigzp)*(zn-z(i)))/(z(i+1)-z(i));
    gs=(u(iigzm)*(z(i)-zs)+u(iig)*(zs-z(i-1)))/(z(i)-z(i-1));

    dgn=(u(iigzp)-u(iig))/(z(i+1)-z(i));
    dgs=(u(iig)-u(iigzm))/(z(i)-z(i-1));

    g(iig)=2*u(iif)*u(iig)+(gn-gs)/(zn-zs)*u(iih)-(dgn-dgs)/(zn-zs);
    g(iig)=g(iig)+2*u(iif);
    
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;2*u(iif)];
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;((z(i+1)-zn))/(z(i+1)-z(i))/(zn-zs)*u(iih)];
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;-((zs-z(i-1)))/(z(i)-z(i-1))/(zn-zs)*u(iih)];
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;(-1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;(1)/(z(i)-z(i-1))*(1)/(zn-zs)];

    ii=[ii;iig]; jj=[jj;iigzp]; vv=[vv;((zn-z(i)))/(z(i+1)-z(i))/(zn-zs)*u(iih)];
    ii=[ii;iig]; jj=[jj;iigzp]; vv=[vv;(1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
    
    ii=[ii;iig]; jj=[jj;iigzm]; vv=[vv;-((z(i)-zs))/(z(i)-z(i-1))/(zn-zs)*u(iih)];
    ii=[ii;iig]; jj=[jj;iigzm]; vv=[vv;(1)/(z(i)-z(i-1))*-(1)/(zn-zs)];

    ii=[ii;iig]; jj=[jj;iif]; vv=[vv;2*u(iig)];
    ii=[ii;iig]; jj=[jj;iih]; vv=[vv;(gn-gs)/(zn-zs)];

    ii=[ii;iig]; jj=[jj;iif]; vv=[vv;2];
    
    end
% 
% 
% %     %zmom WW'+P"=W''
%     hn=(u(iih)*(z(i+1)-zn)+u(iihzp)*(zn-z(i)))/(z(i+1)-z(i));
%     hs=(u(iihzm)*(z(i)-zs)+u(iih)*(zs-z(i-1)))/(z(i)-z(i-1));
% 
%     dhn=(u(iihzp)-u(iih))/(z(i+1)-z(i));
%     dhs=(u(iih)-u(iihzm))/(z(i)-z(i-1));
% 
%     g(iih)=(u(iipzp)-u(iip))/(zn-zs)+(hn-hs)/(zn-zs)*u(iih)-(dhn-dhs)/(zn-zs);
%     
%     
%     ii=[ii;iih]; jj=[jj;iipzp]; vv=[vv;1/(zn-zs)];
%     ii=[ii;iih]; jj=[jj;iip]; vv=[vv;-1/(zn-zs)];
%     
%     ii=[ii;iih]; jj=[jj;iih]; vv=[vv;(hn-hs)/(zn-zs)];
%     ii=[ii;iih]; jj=[jj;iih]; vv=[vv;((z(i+1)-zn))/(z(i+1)-z(i))/(zn-zs)*u(iih)];
%     ii=[ii;iih]; jj=[jj;iih]; vv=[vv;-((zs-z(i-1)))/(z(i)-z(i-1))/(zn-zs)*u(iih)];
%     ii=[ii;iih]; jj=[jj;iih]; vv=[vv;(-1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
%     ii=[ii;iih]; jj=[jj;iih]; vv=[vv;(1)/(z(i)-z(i-1))*(1)/(zn-zs)];
% 
%     ii=[ii;iih]; jj=[jj;iihzp]; vv=[vv;((zn-z(i)))/(z(i+1)-z(i))/(zn-zs)*u(iih)];
%     ii=[ii;iih]; jj=[jj;iihzp]; vv=[vv;(1)/(z(i+1)-z(i))*-(1)/(zn-zs)];
%     
%     ii=[ii;iih]; jj=[jj;iihzm]; vv=[vv;-((z(i)-zs))/(z(i)-z(i-1))/(zn-zs)*u(iih)];
%     ii=[ii;iih]; jj=[jj;iihzm]; vv=[vv;(1)/(z(i)-z(i-1))*-(1)/(zn-zs)];
% 
%   
% end
% 
% %bc
% i=1; iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4; 
% g(iip)=u(iip); g(iif)=u(iif); g(iig)=u(iig); g(iih)=u(iih);
% 
% ii=[ii;iip]; jj=[jj;iip]; vv=[vv;1];
% ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1];
% ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1];
% ii=[ii;iih]; jj=[jj;iih]; vv=[vv;1];
%   
% 
% i=np; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4; 
% g(iif)=u(iif); g(iig)=u(iig)+1; g(iih)=u(5);
% 
% ii=[ii;iih]; jj=[jj;5]; vv=[vv;1];
% ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1];
% ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1];


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