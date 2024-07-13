   %vonK code 
  %p at wall for the p bc which is at wall
%          ----- dummy p
   %       | . | ghost u,v,w
   % z^    ----- p,      | continuity
   %  |    | . | u,v,w | x,y,z mom
   %  |    -----
   %  |    | . |
   %  |    -----
   %  |    | . |
   %  |-----------------> 
%          | . | ghost u,v,w
%          -----



function [g,jac]=evalJacRhs(u,z)

% np=length(z);
zw=[2*z(1)-z(2),z(1:end-1),2*z(end-1)-z(end-2)];
zc=(zw(1:end-1)+zw(2:end))/2;
zw=zw(2:end);
% plot(zc,zc,'-x'); hold on; plot(zw,zw,'-x');
% ii=[; jj=[; vv=[;
ii=ones(length(zw)*40,1)*length(zw)*4; jj=ii; vv=ii*0; iii=1;
g=u*0;

%conti
for i=1:length(zw)-1
    iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
    iifzp=iif+4; iihzp=(i-1)*4+4+4;
    
    zs=zc(i); zn=zc(i+1); zp=zw(i); 
    
%     dz=zn-zs;
    fp=(u(iif)*(zn-zp)+u(iifzp)*(zp-zs))/(zn-zs);    
    g(iip)=2*fp+(u(iihzp)-u(iih))/(zn-zs);

   ii(iii)=iip; jj(iii)=iif;    vv(iii)=2*((zn-zp))/(zn-zs);    iii=iii+1;
   ii(iii)=iip; jj(iii)=iifzp;  vv(iii)=2*((zp-zs))/(zn-zs);    iii=iii+1;
   ii(iii)=iip; jj(iii)=iih;    vv(iii)=-1/(zn-zs);             iii=iii+1;
   ii(iii)=iip; jj(iii)=iihzp;  vv(iii)=1/(zn-zs);              iii=iii+1;
end

%xmom ymom zmom
for i=2:length(zc)-1
    iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
    iifzm=iif-4; iigzm=iig-4; iihzm=iih-4; 
    iifzp=iif+4; iigzp=iig+4; iihzp=iih+4; 
    iipzp=iip+4; iipzm=iip-4;
    
    zp=zc(i); zn=zw(i); zs=zw(i-1); zpm=zc(i-1); zpp=zc(i+1);
    

    %xmom
    fn=(u(iif)*(zpp-zn)+u(iifzp)*(zn-zp))/(zpp-zp);
    fs=(u(iifzm)*(zp-zs)+u(iif)*(zs-zpm))/(zp-zpm);

    dfn=(u(iifzp)-u(iif))/(zpp-zp);
    dfs=(u(iif)-u(iifzm))/(zp-zpm);

    g(iif)=u(iif)^2-u(iig)^2+(fn-fs)/(zn-zs)*u(iih)-(dfn-dfs)/(zn-zs);
%     % g(iig)=g(iig)-u(iig)^2-2*u(iig)-1;
    g(iif)=g(iif)-2*u(iig)-1;
%    
%     % 
    ii(iii)=iif; jj(iii)=iif; vv(iii)=2*u(iif);                                 iii=iii+1;
    ii(iii)=iif; jj(iii)=iif; vv(iii)=((zpp-zn))/(zpp-zp)/(zn-zs)*u(iih);       iii=iii+1;
    ii(iii)=iif; jj(iii)=iif; vv(iii)=-((zs-zpm))/(zp-zpm)/(zn-zs)*u(iih);      iii=iii+1;
    ii(iii)=iif; jj(iii)=iif; vv(iii)=(-1)/(zpp-zp)*-(1)/(zn-zs);               iii=iii+1;
    ii(iii)=iif; jj(iii)=iif; vv(iii)=(1)/(zp-zpm)*(1)/(zn-zs);                 iii=iii+1;

    ii(iii)=iif; jj(iii)=iifzp; vv(iii)=((zn-zp))/(zpp-zp)/(zn-zs)*u(iih);      iii=iii+1;
    ii(iii)=iif; jj(iii)=iifzp; vv(iii)=(1)/(zpp-zp)*-(1)/(zn-zs);              iii=iii+1;
    
    ii(iii)=iif; jj(iii)=iifzm; vv(iii)=-((zp-zs))/(zp-zpm)/(zn-zs)*u(iih);     iii=iii+1;
    ii(iii)=iif; jj(iii)=iifzm; vv(iii)=(1)/(zp-zpm)*-(1)/(zn-zs);              iii=iii+1;
% 
    ii(iii)=iif; jj(iii)=iig; vv(iii)=-2*u(iig);                                iii=iii+1;
    ii(iii)=iif; jj(iii)=iih; vv(iii)=(fn-fs)/(zn-zs);                          iii=iii+1;
% 
    ii(iii)=iif; jj(iii)=iig; vv(iii)=-2;                                       iii=iii+1;
% 
%     %ymom
    gn=(u(iig)*(zpp-zn)+u(iigzp)*(zn-zp))/(zpp-zp);
    gs=(u(iigzm)*(zp-zs)+u(iig)*(zs-zpm))/(zp-zpm);

    dgn=(u(iigzp)-u(iig))/(zpp-zp);
    dgs=(u(iig)-u(iigzm))/(zp-zpm);

    g(iig)=2*u(iif)*u(iig)+(gn-gs)/(zn-zs)*u(iih)-(dgn-dgs)/(zn-zs);
    g(iig)=g(iig)+2*u(iif);
    
    ii(iii)=iig; jj(iii)=iig; vv(iii)=2*u(iif); iii=iii+1;
    ii(iii)=iig; jj(iii)=iig; vv(iii)=((zpp-zn))/(zpp-zp)/(zn-zs)*u(iih);iii=iii+1;
    ii(iii)=iig; jj(iii)=iig; vv(iii)=-((zs-zpm))/(zp-zpm)/(zn-zs)*u(iih);iii=iii+1;
    ii(iii)=iig; jj(iii)=iig; vv(iii)=(-1)/(zpp-zp)*-(1)/(zn-zs);iii=iii+1;
    ii(iii)=iig; jj(iii)=iig; vv(iii)=(1)/(zp-zpm)*(1)/(zn-zs);iii=iii+1;

    ii(iii)=iig; jj(iii)=iigzp; vv(iii)=((zn-zp))/(zpp-zp)/(zn-zs)*u(iih);iii=iii+1;
    ii(iii)=iig; jj(iii)=iigzp; vv(iii)=(1)/(zpp-zp)*-(1)/(zn-zs);iii=iii+1;
    
    ii(iii)=iig; jj(iii)=iigzm; vv(iii)=-((zp-zs))/(zp-zpm)/(zn-zs)*u(iih);iii=iii+1;
    ii(iii)=iig; jj(iii)=iigzm; vv(iii)=(1)/(zp-zpm)*-(1)/(zn-zs);iii=iii+1;

    ii(iii)=iig; jj(iii)=iif; vv(iii)=2*u(iig);iii=iii+1;
    ii(iii)=iig; jj(iii)=iih; vv(iii)=(gn-gs)/(zn-zs);iii=iii+1;

    ii(iii)=iig; jj(iii)=iif; vv(iii)=2;iii=iii+1;
% 
% 
% %     %zmom WW'+P"=W''
    hn=(u(iih)*(zpp-zn)+u(iihzp)*(zn-zp))/(zpp-zp);
    hs=(u(iihzm)*(zp-zs)+u(iih)*(zs-zpm))/(zp-zpm);

    dhn=(u(iihzp)-u(iih))/(zpp-zp);
    dhs=(u(iih)-u(iihzm))/(zp-zpm);

    g(iih)=(u(iip)-u(iipzm))/(zn-zs)+(hn-hs)/(zn-zs)*u(iih)-(dhn-dhs)/(zn-zs);
    
    
    ii(iii)=iih; jj(iii)=iip; vv(iii)=1/(zn-zs);iii=iii+1;
    ii(iii)=iih; jj(iii)=iipzm; vv(iii)=-1/(zn-zs);iii=iii+1;
    
    ii(iii)=iih; jj(iii)=iih; vv(iii)=(hn-hs)/(zn-zs);iii=iii+1;
    ii(iii)=iih; jj(iii)=iih; vv(iii)=((zpp-zn))/(zpp-zp)/(zn-zs)*u(iih);iii=iii+1;
    ii(iii)=iih; jj(iii)=iih; vv(iii)=-((zs-zpm))/(zp-zpm)/(zn-zs)*u(iih);iii=iii+1;
    ii(iii)=iih; jj(iii)=iih; vv(iii)=(-1)/(zpp-zp)*-(1)/(zn-zs);iii=iii+1;
    ii(iii)=iih; jj(iii)=iih; vv(iii)=(1)/(zp-zpm)*(1)/(zn-zs);iii=iii+1;

    ii(iii)=iih; jj(iii)=iihzp; vv(iii)=((zn-zp))/(zpp-zp)/(zn-zs)*u(iih);iii=iii+1;
    ii(iii)=iih; jj(iii)=iihzp; vv(iii)=(1)/(zpp-zp)*-(1)/(zn-zs);iii=iii+1;
    
    ii(iii)=iih; jj(iii)=iihzm; vv(iii)=-((zp-zs))/(zp-zpm)/(zn-zs)*u(iih);iii=iii+1;
    ii(iii)=iih; jj(iii)=iihzm; vv(iii)=(1)/(zp-zpm)*-(1)/(zn-zs);iii=iii+1;

  
end

%bc
i=1; iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4; 
iifzp=iif+4;  iigzp=iig+4;  iihzp=iih+4; 
% g(iip)=u(iip);
g(iif)=(u(iif)+u(iifzp))/2; g(iig)=(u(iig)+u(iigzp))/2; g(iih)=(u(iih)+u(iihzp))/2;

% ii(iii)=iip; jj(iii)=iip; vv(iii)=1;
ii(iii)=iif; jj(iii)=iif; vv(iii)=1/2;iii=iii+1;
ii(iii)=iif; jj(iii)=iifzp; vv(iii)=1/2;iii=iii+1;
ii(iii)=iig; jj(iii)=iig; vv(iii)=1/2;iii=iii+1;
ii(iii)=iig; jj(iii)=iigzp; vv(iii)=1/2;iii=iii+1;
ii(iii)=iih; jj(iii)=iih; vv(iii)=1/2;iii=iii+1;
ii(iii)=iih; jj(iii)=iihzp; vv(iii)=1/2;iii=iii+1;
%   
% 
i=length(zw); iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
iifzm=iif-4; iigzm=iig-4;
g(iif)=(u(iif)+u(iifzm))/2; g(iig)=(u(iig)+u(iigzm))/2+1; g(iih)=u(1); g(iip)=u(iip)-u(iip-4);

ii(iii)=iih; jj(iii)=1; vv(iii)=1;iii=iii+1;
ii(iii)=iif; jj(iii)=iif; vv(iii)=1/2;iii=iii+1;
ii(iii)=iif; jj(iii)=iifzm; vv(iii)=1/2;iii=iii+1;
ii(iii)=iig; jj(iii)=iig; vv(iii)=1/2;iii=iii+1;
ii(iii)=iig; jj(iii)=iigzm; vv(iii)=1/2;iii=iii+1;
ii(iii)=iip; jj(iii)=iip; vv(iii)=1;iii=iii+1;
ii(iii)=iip; jj(iii)=iip-4; vv(iii)=-1;iii=iii+1;


jac=sparse(ii,jj,vv);
end

% function dydx = bvpfcn(x,y) % equation to solve
% dydx = zeros(5,1);
% dydx = [y(4)
%        y(5)
%        -2*y(1)
%        y(1)^2-y(2)^2+y(4)*y(3)
%        2*y(1)*y(2)+y(5)*y(3);
% end
% %--------------------------------
% function res = bcfcn(ya,yb) % boundary conditions
% res = [ya(1)
%     ya(2)-1
%     ya(3)
%     yb(1)
%     yb(2) ;
% end