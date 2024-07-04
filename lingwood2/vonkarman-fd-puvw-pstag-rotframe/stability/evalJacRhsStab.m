%vonK code
%uvw at cell wall, p at cell center, mom eqs written at wall,
%conti at center

% z^    ----- u,v,w | x,y,z mom
%  |    | . | p,      | continuity
%  |    -----
%  |    | . |
%  |    -----
%  |    | . |
%  |----------------->
%          | . | dummy p
%          -----



function [g,jac0,jac1]=evalJacRhsStab(u,z,U,omega,bbar,R,alpha)

np=length(z);
ii=[np*4]; jj=[np*4]; vv=[0];
ii1=[np*4]; jj1=[np*4]; vv1=[0];
g=u*0;

%conti
for i=2:np
    iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
    iifzm=iif-4; iigzm=iig-4; iihzm=(i-1)*4+4-4;
    dz=z(i)-z(i-1);
    
    g(iip)=1i*bbar*(u(iig)+u(iigzm))/2+(u(iih)-u(iihzm))/dz;
    
    ii=[ii;iip]; jj=[jj;iih]; vv=[vv;1/dz];
    ii=[ii;iip]; jj=[jj;iihzm]; vv=[vv;-1/dz];
    
    ii=[ii;iip]; jj=[jj;iig]; vv=[vv;1i*bbar/2];
    ii=[ii;iip]; jj=[jj;iigzm]; vv=[vv;1i*bbar/2];
    
    
    g(iip)=g(iip)+1i*alpha*(u(iif)+u(iifzm))/2;
    
    ii1=[ii1;iip]; jj1=[jj1;iif]; vv1=[vv1;1i/2];
    ii1=[ii1;iip]; jj1=[jj1;iifzm]; vv1=[vv1;1i/2];
    
end

%xmom ymom zmom
for i=2:np-1
    iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
    iifzm=iif-4; iigzm=iig-4; iihzm=iih-4;
    iifzp=iif+4; iigzp=iig+4; iihzp=iih+4;
    iipzp=iip+4;
    
    zp=z(i); zn=(z(i)+z(i+1))/2; zs=(z(i)+z(i-1))/2;
    
    
    %xmom
    %     fn=(u(iif)*(z(i+1)-zn)+u(iifzp)*(zn-z(i)))/(z(i+1)-z(i));
    %     fs=(u(iifzm)*(z(i)-zs)+u(iif)*(zs-z(i-1)))/(z(i)-z(i-1));
    Fn=(U(iif)*(z(i+1)-zn)+U(iifzp)*(zn-z(i)))/(z(i+1)-z(i));
    Fs=(U(iifzm)*(z(i)-zs)+U(iif)*(zs-z(i-1)))/(z(i)-z(i-1));
    pp=(u(iip)*(zn-z(i))+u(iipzp)*(z(i)-zs))/(zn-zs);
    
    %     dfn=(u(iifzp)-u(iif))/(z(i+1)-z(i));
    %     dfs=(u(iif)-u(iifzm))/(z(i)-z(i-1));
    
    g(iif)=u(iif)*(1i*bbar*U(iig)-1i*omega)+u(iih)*(Fn-Fs)/(zn-zs);
    g(iif)=g(iif)+u(iif)*1i*alpha*U(iif)+1i*alpha*pp;
    
    %
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;(1i*bbar*U(iig)-1i*omega)];
    ii=[ii;iif]; jj=[jj;iih]; vv=[vv;(Fn-Fs)/(zn-zs)];
    
    ii1=[ii1;iif]; jj1=[jj1;iif]; vv1=[vv1;1i*U(iif)];
    ii1=[ii1;iif]; jj1=[jj1;iip]; vv1=[vv1;1i*((zn-z(i)))/(zn-zs)];
    ii1=[ii1;iif]; jj1=[jj1;iipzp]; vv1=[vv1;1i*((z(i)-zs))/(zn-zs)];
    
    %y mom
    Gn=(U(iig)*(z(i+1)-zn)+U(iigzp)*(zn-z(i)))/(z(i+1)-z(i));
    Gs=(U(iigzm)*(z(i)-zs)+U(iig)*(zs-z(i-1)))/(z(i)-z(i-1));
    pp=(u(iip)*(zn-z(i))+u(iipzp)*(z(i)-zs))/(zn-zs);
    
    %     dfn=(u(iifzp)-u(iif))/(z(i+1)-z(i));
    %     dfs=(u(iif)-u(iifzm))/(z(i)-z(i-1));
    
    g(iig)=u(iig)*(1i*bbar*U(iig)-1i*omega)+u(iih)*(Gn-Gs)/(zn-zs)+1i*bbar*pp;
    g(iig)=g(iig)+u(iig)*1i*alpha*U(iif);
    
    %
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;(1i*bbar*U(iig)-1i*omega)];
    ii=[ii;iig]; jj=[jj;iih]; vv=[vv;(Gn-Gs)/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iip]; vv=[vv;bbar*1i*((zn-z(i)))/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iipzp]; vv=[vv;bbar*1i*((z(i)-zs))/(zn-zs)];
    
    
    ii1=[ii1;iig]; jj1=[jj1;iig]; vv1=[vv1;1i*U(iif)];
    
    %zmom
    g(iih)=u(iih)*(1i*bbar*U(iig)-1i*omega)+(u(iipzp)-u(iip))/(zn-zs);
    g(iih)=g(iih)+u(iih)*1i*alpha*U(iif);
    
    %
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;(1i*bbar*U(iig)-1i*omega)];
    ii=[ii;iih]; jj=[jj;iip]; vv=[vv;-1/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iipzp]; vv=[vv;1/(zn-zs)];
    
    ii1=[ii1;iih]; jj1=[jj1;iih]; vv1=[vv1;1i*U(iif)];
    
end

%bc
i=1; iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
g(iip)=u(iip); g(iif)=u(iif); g(iig)=u(iig); g(iih)=u(iih);

ii=[ii;iip]; jj=[jj;iip]; vv=[vv;1];
ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1];
ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1];
ii=[ii;iih]; jj=[jj;iih]; vv=[vv;1];


i=np; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
g(iif)=u(iif); g(iig)=u(iig); g(iih)=u(5);

ii=[ii;iih]; jj=[jj;5]; vv=[vv;1];
ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1];
ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1];


jac0=sparse(ii,jj,vv);
jac1=sparse(ii1,jj1,vv1);
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