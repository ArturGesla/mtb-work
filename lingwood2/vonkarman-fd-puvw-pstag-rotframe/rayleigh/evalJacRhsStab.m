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



function [g,jac0,jac1,jac2,jac3]=evalJacRhsStab(u,z,U,omega,bbar,R,alpha)

np=length(z);
ii=[np]; jj=[np]; vv=[0];
ii1=[np]; jj1=[np]; vv1=[0];
ii2=[np]; jj2=[np]; vv2=[0];
ii3=[np]; jj3=[np]; vv3=[0];

g=zeros(np,1);

%rayleigh
for i=2:np-1
%     iip=(i-1)*4+1; 
    iiF=(i-1)*4+2; iiG=(i-1)*4+3; 
    iih=(i-1)*1+1;
    iiFzm=iiF-4; iiGzm=iiG-4; 
    iihzm=iih-1;
    iiFzp=iiF+4; iiGzp=iiG+4;
    iihzp=iih+1;
%     iipzp=iip+4;
    
    zp=z(i); zn=(z(i)+z(i+1))/2; zs=(z(i)+z(i-1))/2;
    
    
    %xmom
    %     fn=(u(iif)*(z(i+1)-zn)+u(iifzp)*(zn-z(i)))/(z(i+1)-z(i));
    %     fs=(u(iifzm)*(z(i)-zs)+u(iif)*(zs-z(i-1)))/(z(i)-z(i-1));
%     Fn=(U(iif)*(z(i+1)-zn)+U(iifzp)*(zn-z(i)))/(z(i+1)-z(i));
%     Fs=(U(iifzm)*(z(i)-zs)+U(iif)*(zs-z(i-1)))/(z(i)-z(i-1));
%     pp=(u(iip)*(zn-z(i))+u(iipzp)*(z(i)-zs))/(zn-zs);
    
        dFn=(U(iiFzp)-U(iiF))/(z(i+1)-z(i));
        dFs=(U(iiF)-U(iiFzm))/(z(i)-z(i-1));
        dGn=(U(iiGzp)-U(iiG))/(z(i+1)-z(i));
        dGs=(U(iiG)-U(iiGzm))/(z(i)-z(i-1));
        dhn=(u(iihzp)-u(iih))/(z(i+1)-z(i));
        dhs=(u(iih)-u(iihzm))/(z(i)-z(i-1));
    
    
    g(iih)=g(iih)+alpha*U(iiF)*(dhn-dhs)/(zn-zs);
    ii1=[ii1;iih]; jj1=[jj1;iih]; vv1=[vv1;-U(iiF)/(z(i+1)-z(i))/(zn-zs)];
    ii1=[ii1;iih]; jj1=[jj1;iih]; vv1=[vv1;-U(iiF)/(z(i)-z(i-1))/(zn-zs)];
    ii1=[ii1;iih]; jj1=[jj1;iihzp]; vv1=[vv1;U(iiF)/(z(i+1)-z(i))/(zn-zs)];
    ii1=[ii1;iih]; jj1=[jj1;iihzm]; vv1=[vv1;U(iiF)/(z(i)-z(i-1))/(zn-zs)];
    
    g(iih)=g(iih)-alpha.^3*U(iiF)*u(iih);
    ii3=[ii3;iih]; jj3=[jj3;iih]; vv3=[vv3;-1*U(iiF)];
    
    g(iih)=g(iih)-alpha*bbar.^2*U(iiF)*u(iih);
    ii1=[ii1;iih]; jj1=[jj1;iih]; vv1=[vv1;-bbar.^2*U(iiF)];
    
    g(iih)=g(iih)+(bbar*U(iiG)-omega)*(dhn-dhs)/(zn-zs);
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;-(bbar*U(iiG)-omega)/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;-(bbar*U(iiG)-omega)/(z(i)-z(i-1))/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iihzp]; vv=[vv;(bbar*U(iiG)-omega)/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iihzm]; vv=[vv;(bbar*U(iiG)-omega)/(z(i)-z(i-1))/(zn-zs)];
    
    g(iih)=g(iih)-(bbar*U(iiG)-omega)*alpha.^2*u(iih);
    ii2=[ii2;iih]; jj2=[jj2;iih]; vv2=[vv2;-(bbar*U(iiG)-omega)];
    
    g(iih)=g(iih)-(bbar*U(iiG)-omega)*bbar.^2*u(iih);
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;-(bbar*U(iiG)-omega)*bbar.^2];
    
    g(iih)=g(iih)-alpha*(dFn-dFs)/(zn-zs)*u(iih);
    ii1=[ii1;iih]; jj1=[jj1;iih]; vv1=[vv1;-(dFn-dFs)/(zn-zs)];
    
    g(iih)=g(iih)-bbar*(dGn-dGs)/(zn-zs)*u(iih);
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;-bbar*(dGn-dGs)/(zn-zs)];
    
    
    
end

%bc
i=1; 
iih=(i-1)*1+1;    iihzm=iih-1;       iihzp=iih+1;

g(iih)=u(iih); 
ii=[ii;iih]; jj=[jj;iih]; vv=[vv;1];


% %dirichlet at inf
% i=np; 
% iih=(i-1)*1+1;    iihzm=iih-1;       iihzp=iih+1;
% g(iih)=u(iih); 
% ii=[ii;iih]; jj=[jj;iih]; vv=[vv;1];

%neumann at inf
i=np-1; 
iih=(i-1)*1+1;    iihzm=iih-1;       iihzp=iih+1;
hn=(u(iih)*(z(i+1)-zn)+u(iihzp)*(zn-z(i)))/(z(i+1)-z(i));
hs=(u(iihzm)*(z(i)-zs)+u(iih)*(zs-z(i-1)))/(z(i)-z(i-1));
zp=z(i); zn=(z(i)+z(i+1))/2; zs=(z(i)+z(i-1))/2;

g(iihzp)=(hn-hs)/(zn-zs); 
ii=[ii;iihzp]; jj=[jj;iih]; vv=[vv;(z(i+1)-zn)/(z(i+1)-z(i))/(zn-zs)];
ii=[ii;iihzp]; jj=[jj;iih]; vv=[vv;-(zs-z(i-1))/(z(i)-z(i-1))/(zn-zs)];
ii=[ii;iihzp]; jj=[jj;iihzp]; vv=[vv;(zn-z(i))/(z(i+1)-z(i))/(zn-zs)];
ii=[ii;iihzp]; jj=[jj;iihzm]; vv=[vv;-(z(i)-zs)/(z(i)-z(i-1))/(zn-zs)];

% 
% 
% i=np; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
% g(iif)=u(iif); g(iig)=u(iig); g(iih)=u(5);
% 
% ii=[ii;iih]; jj=[jj;5]; vv=[vv;1];
% ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1];
% ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1];


jac0=sparse(ii,jj,vv);
jac1=sparse(ii1,jj1,vv1);
jac2=sparse(ii2,jj2,vv2);
jac3=sparse(ii3,jj3,vv3);
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