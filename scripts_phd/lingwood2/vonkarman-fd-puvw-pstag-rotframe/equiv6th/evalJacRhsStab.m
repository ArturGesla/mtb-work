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



function [g,jac0,jac1,jac2]=evalJacRhsStab(u,z,U,omega,bbar,R,alpha)

np=length(z);
ii=[np*4]; jj=[np*4]; vv=[0];
ii1=[np*4]; jj1=[np*4]; vv1=[0];
ii2=[np*4]; jj2=[np*4]; vv2=[0];
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
    
    %1/R parts 
    g(iip)=g(iip)+(u(iif)+u(iifzm))/2/R;
    
    ii=[ii;iip]; jj=[jj;iif]; vv=[vv;1/2/R];
    ii=[ii;iip]; jj=[jj;iifzm]; vv=[vv;1/2/R];
    
end

%xmom ymom zmom
for i=2:np-1
    iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
    iifzm=iif-4; iigzm=iig-4; iihzm=iih-4;
    iifzp=iif+4; iigzp=iig+4; iihzp=iih+4;
    iipzp=iip+4;
    
    zp=z(i); zn=(z(i)+z(i+1))/2; zs=(z(i)+z(i-1))/2;
    
    
    %xmom
        fn=(u(iif)*(z(i+1)-zn)+u(iifzp)*(zn-z(i)))/(z(i+1)-z(i));
        fs=(u(iifzm)*(z(i)-zs)+u(iif)*(zs-z(i-1)))/(z(i)-z(i-1));
    Fn=(U(iif)*(z(i+1)-zn)+U(iifzp)*(zn-z(i)))/(z(i+1)-z(i));
    Fs=(U(iifzm)*(z(i)-zs)+U(iif)*(zs-z(i-1)))/(z(i)-z(i-1));
    pp=(u(iip)*(zn-z(i))+u(iipzp)*(z(i)-zs))/(zn-zs);
    
        dfn=(u(iifzp)-u(iif))/(z(i+1)-z(i));
        dfs=(u(iif)-u(iifzm))/(z(i)-z(i-1));
    
    g(iif)=u(iif)*(1i*bbar*U(iig)-1i*omega)+u(iih)*(Fn-Fs)/(zn-zs);
    g(iif)=g(iif)+u(iif)*1i*alpha*U(iif)+1i*alpha*pp;
    
    %
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;(1i*bbar*U(iig)-1i*omega)];
    ii=[ii;iif]; jj=[jj;iih]; vv=[vv;(Fn-Fs)/(zn-zs)];
    
    ii1=[ii1;iif]; jj1=[jj1;iif]; vv1=[vv1;1i*U(iif)];
    ii1=[ii1;iif]; jj1=[jj1;iip]; vv1=[vv1;1i*((zn-z(i)))/(zn-zs)];
    ii1=[ii1;iif]; jj1=[jj1;iipzp]; vv1=[vv1;1i*((z(i)-zs))/(zn-zs)];
    
    %1/R parts
    g(iif)=g(iif)-1/R*(dfn-dfs)/(zn-zs);
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1/R*1/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1/R*1/(z(i)-z(i-1))/(zn-zs)];
    ii=[ii;iif]; jj=[jj;iifzp]; vv=[vv;1/R*-1/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iif]; jj=[jj;iifzm]; vv=[vv;1/R*-1/(z(i)-z(i-1))/(zn-zs)];
    
    g(iif)=g(iif)+1/R*u(iif)*U(iif)+1/R*U(iih)*(fn-fs)/(zn-zs)-2/R*U(iig)*u(iig)-2/R*u(iig)+u(iif)*bbar^2/R;
    
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1/R*U(iif)];
    
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1/R*U(iih)*((z(i+1)-zn))/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iif]; jj=[jj;iifzp]; vv=[vv;1/R*U(iih)*((zn-z(i)))/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iif]; jj=[jj;iifzm]; vv=[vv;-1/R*U(iih)*((z(i)-zs))/(z(i)-z(i-1))/(zn-zs)];
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;-1/R*U(iih)*((zs-z(i-1)))/(z(i)-z(i-1))/(zn-zs)];
    
    ii=[ii;iif]; jj=[jj;iig]; vv=[vv;-2/R*U(iig)];
    ii=[ii;iif]; jj=[jj;iig]; vv=[vv;-2/R];
    ii=[ii;iif]; jj=[jj;iif]; vv=[vv;bbar^2/R];
    
    
    g(iif)=g(iif)+u(iif)*alpha^2/R;
    ii2=[ii2;iif]; jj2=[jj2;iif]; vv2=[vv2;1/R];

    
    
    %y mom
    gn=(u(iig)*(z(i+1)-zn)+u(iigzp)*(zn-z(i)))/(z(i+1)-z(i));
    gs=(u(iigzm)*(z(i)-zs)+u(iig)*(zs-z(i-1)))/(z(i)-z(i-1));
    Gn=(U(iig)*(z(i+1)-zn)+U(iigzp)*(zn-z(i)))/(z(i+1)-z(i));
    Gs=(U(iigzm)*(z(i)-zs)+U(iig)*(zs-z(i-1)))/(z(i)-z(i-1));
    pp=(u(iip)*(zn-z(i))+u(iipzp)*(z(i)-zs))/(zn-zs);
    
        dgn=(u(iigzp)-u(iig))/(z(i+1)-z(i));
        dgs=(u(iig)-u(iigzm))/(z(i)-z(i-1));
    
    g(iig)=u(iig)*(1i*bbar*U(iig)-1i*omega)+u(iih)*(Gn-Gs)/(zn-zs)+1i*bbar*pp;
    g(iig)=g(iig)+u(iig)*1i*alpha*U(iif);
    
    %
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;(1i*bbar*U(iig)-1i*omega)];
    ii=[ii;iig]; jj=[jj;iih]; vv=[vv;(Gn-Gs)/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iip]; vv=[vv;bbar*1i*((zn-z(i)))/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iipzp]; vv=[vv;bbar*1i*((z(i)-zs))/(zn-zs)];
    
    ii1=[ii1;iig]; jj1=[jj1;iig]; vv1=[vv1;1i*U(iif)];
    
    %1/R terms
    g(iig)=g(iig)-1/R*(dgn-dgs)/(zn-zs);
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1/R*1/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1/R*1/(z(i)-z(i-1))/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iigzp]; vv=[vv;1/R*-1/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iigzm]; vv=[vv;1/R*-1/(z(i)-z(i-1))/(zn-zs)];
    
    g(iig)=g(iig)+2/R*u(iif)*U(iig)+1/R*U(iif)*u(iig)+1/R*U(iih)*(gn-gs)/(zn-zs)+2/R*u(iif)+1/R*bbar^2*u(iig);
    
    ii=[ii;iig]; jj=[jj;iif]; vv=[vv;2/R*U(iig)];
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1/R*U(iif)];
    
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1/R*U(iih)*((z(i+1)-zn))/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iigzp]; vv=[vv;1/R*U(iih)*((zn-z(i)))/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;-1/R*U(iih)*((zs-z(i-1)))/(z(i)-z(i-1))/(zn-zs)];
    ii=[ii;iig]; jj=[jj;iigzm]; vv=[vv;-1/R*U(iih)*((z(i)-zs))/(z(i)-z(i-1))/(zn-zs)];
    
    ii=[ii;iig]; jj=[jj;iif]; vv=[vv;2/R];
    ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1/R*bbar^2];
    
    

    g(iig)=g(iig)+1/R*alpha^2*u(iig);
    ii2=[ii2;iig]; jj2=[jj2;iig]; vv2=[vv2;1/R];
    
    
    
    %zmom
    Hn=(U(iih)*(z(i+1)-zn)+U(iihzp)*(zn-z(i)))/(z(i+1)-z(i));
    Hs=(U(iihzm)*(z(i)-zs)+U(iih)*(zs-z(i-1)))/(z(i)-z(i-1));
    
    hn=(u(iih)*(z(i+1)-zn)+u(iihzp)*(zn-z(i)))/(z(i+1)-z(i));
    hs=(u(iihzm)*(z(i)-zs)+u(iih)*(zs-z(i-1)))/(z(i)-z(i-1));
    
    dhn=(u(iihzp)-u(iih))/(z(i+1)-z(i));
    dhs=(u(iih)-u(iihzm))/(z(i)-z(i-1));
    
    g(iih)=u(iih)*(1i*bbar*U(iig)-1i*omega)+(u(iipzp)-u(iip))/(zn-zs);
    
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;(1i*bbar*U(iig)-1i*omega)];
    ii=[ii;iih]; jj=[jj;iip]; vv=[vv;-1/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iipzp]; vv=[vv;1/(zn-zs)];
    
    g(iih)=g(iih)+u(iih)*1i*alpha*U(iif);
    ii1=[ii1;iih]; jj1=[jj1;iih]; vv1=[vv1;1i*U(iif)];
    
    %1/R parts
    g(iih)=g(iih)+u(iih)*alpha^2/R;
    ii2=[ii2;iih]; jj2=[jj2;iih]; vv2=[vv2; 1/R];
    
    g(iih)=g(iih)+u(iih)*(bbar.^2/R+1/R*(Hn-Hs)/(zn-zs));
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;(bbar.^2/R+1/R*(Hn-Hs)/(zn-zs))];
    
    g(iih)=g(iih)-1/R*(dhn-dhs)/(zn-zs)+1/R*U(iih)*(hn-hs)/(zn-zs);
    
    ii=[ii;iih]; jj=[jj;iihzp]; vv=[vv;-1/R*1/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;-1/R*-1/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;-1*-1/R*1/(z(i)-z(i-1))/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iihzm]; vv=[vv;-1*-1/R*-1/(z(i)-z(i-1))/(zn-zs)];

    ii=[ii;iih]; jj=[jj;iihzp]; vv=[vv;1/R*U(iih)*((zn-z(i)))/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;1/R*U(iih)*((z(i+1)-zn))/(z(i+1)-z(i))/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iihzm]; vv=[vv;-1/R*U(iih)*((z(i)-zs))/(z(i)-z(i-1))/(zn-zs)];
    ii=[ii;iih]; jj=[jj;iih]; vv=[vv;-1/R*U(iih)*((zs-z(i-1)))/(z(i)-z(i-1))/(zn-zs)];
    
end

%bc
i=1; iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
g(iip)=u(iip); g(iif)=u(iif); g(iig)=u(iig); g(iih)=u(iih);

ii=[ii;iip]; jj=[jj;iip]; vv=[vv;1];
ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1];
ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1];
ii=[ii;iih]; jj=[jj;iih]; vv=[vv;1];


i=np; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
g(iif)=u(iif); g(iig)=u(iig); 

% g(iih)=u(5); ii=[ii;iih]; jj=[jj;5]; vv=[vv;1]; %this generates rank-1 def jacobian
g(iih)=u(iih); ii=[ii;iih]; jj=[jj;iih]; vv=[vv;1]; %artificial contion but helped a bit, maybe we should fix p to 0 at 0,inf also
ii=[ii;iif]; jj=[jj;iif]; vv=[vv;1];
ii=[ii;iig]; jj=[jj;iig]; vv=[vv;1];


jac0=sparse(ii,jj,vv);
jac1=sparse(ii1,jj1,vv1);
jac2=sparse(ii2,jj2,vv2);
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