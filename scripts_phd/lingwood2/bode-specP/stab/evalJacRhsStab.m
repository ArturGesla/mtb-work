%vonK code

% z^    ----- p,      | continuity
%  |    | . | u,v,w | x,y,z mom
%  |    -----
%  |    | . |
%  |    -----
%  |    | . |
%  |----------------->
%          | . | 
%          -----



function [g,jac0,jac1,jac2,jacom]=evalJacRhsStab(u,z,U,omega,bbar,R,alpha)

zw=[2*z(1)-z(2),z(1:end-1),2*z(end-1)-z(end-2)];
zc=(zw(1:end-1)+zw(2:end))/2;
zw=zw(2:end);
% plot(zc,zc,'-x'); hold on; plot(zw,zw,'-x');
% ii=[; jj=[; vv=[;
ii=ones(length(zw)*80,1)*length(zw)*4; jj=ii; vv=ii*0; iii=1;
ii1=ones(length(zw)*40/5,1)*length(zw)*4; jj1=ii1; vv1=ii1*0; iii1=1;
ii2=ones(length(zw)*4,1)*length(zw)*4; jj2=ii2; vv2=ii2*0; iii2=1;
iio=ones(length(zw)*4,1)*length(zw)*4; jjo=iio; vvo=iio*0; iiio=1;

g=u*0;

%conti
for i=1:length(zw)-1
%     for i=2:length(zw)-1-1
    iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
    iifzp=iif+4; iigzp=iig+4; iihzp=(i-1)*4+4+4;
    
    
        zs=zc(i); zn=zc(i+1); zp=zw(i); 
        dz=zn-zs;
        
        fp=(u(iif)*(zn-zp)+u(iifzp)*(zp-zs))/(zn-zs);
        gp=(u(iig)*(zn-zp)+u(iigzp)*(zp-zs))/(zn-zs);

    
    
    
    g(iip)=1i*bbar*gp+(u(iihzp)-u(iih))/dz;
    
    ii(iii)=iip; jj(iii)=iihzp; vv(iii)=1/dz; iii=iii+1;
    ii(iii)=iip; jj(iii)=iih; vv(iii)=-1/dz;iii=iii+1;
    
    ii(iii)=iip; jj(iii)=iig; vv(iii)=1i*bbar*((zn-zp))/(zn-zs);iii=iii+1;
    ii(iii)=iip; jj(iii)=iigzp; vv(iii)=1i*bbar*((zp-zs))/(zn-zs);iii=iii+1;
    
    g(iip)=g(iip)+1i*alpha*fp;
    
    ii1(iii1)=iip; jj1(iii1)=iif; vv1(iii1)=1i*((zn-zp))/(zn-zs);iii1=iii1+1;
    ii1(iii1)=iip; jj1(iii1)=iifzp; vv1(iii1)=1i*((zp-zs))/(zn-zs);iii1=iii1+1;
    
    %1/R parts 
    g(iip)=g(iip)+fp/R;
    
    ii(iii)=iip; jj(iii)=iif; vv(iii)=((zn-zp))/(zn-zs)/R;iii=iii+1;
    ii(iii)=iip; jj(iii)=iifzp; vv(iii)=((zp-zs))/(zn-zs)/R;iii=iii+1;
    
end

%xmom ymom zmom
for i=2:length(zc)-1
    iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
    iifzm=iif-4; iigzm=iig-4; iihzm=iih-4;
    iifzp=iif+4; iigzp=iig+4; iihzp=iih+4;
    iipzm=iip-4;
    
    zp=zc(i); zn=zw(i); zs=zw(i-1); zpm=zc(i-1); zpp=zc(i+1);
%     
%     
%     %xmom
        fn=(u(iif)*(zpp-zn)+u(iifzp)*(zn-zp))/(zpp-zp);
        fs=(u(iifzm)*(zp-zs)+u(iif)*(zs-zpm))/(zp-zpm);
    Fn=(U(iif)*(zpp-zn)+U(iifzp)*(zn-zp))/(zpp-zp);
    Fs=(U(iifzm)*(zp-zs)+U(iif)*(zs-zpm))/(zp-zpm);
    pp=(u(iip)+u(iipzm))/2;
    
        dfn=(u(iifzp)-u(iif))/(zpp-zp);
        dfs=(u(iif)-u(iifzm))/(zp-zpm);
    
    g(iif)=u(iif)*(1i*bbar*U(iig)-1i*omega)+u(iih)*(Fn-Fs)/(zn-zs);
    g(iif)=g(iif)+u(iif)*1i*alpha*U(iif)+1i*alpha*pp;
    
    %
    iio(iiio)=iif; jjo(iiio)=iif; vvo(iiio)=-1i; iiio=iiio+1;
    
%     ii(iii)=iif; jj(iii)=iif; vv(iii)=(1i*bbar*U(iig)-1i*omega); iii=iii+1;
    ii(iii)=iif; jj(iii)=iif; vv(iii)=(1i*bbar*U(iig)); iii=iii+1;
    ii(iii)=iif; jj(iii)=iih; vv(iii)=(Fn-Fs)/(zn-zs);iii=iii+1;
    
    ii1(iii1)=iif; jj1(iii1)=iif; vv1(iii1)=1i*U(iif); iii1=iii1+1;
    ii1(iii1)=iif; jj1(iii1)=iip; vv1(iii1)=1i/2;iii1=iii1+1;
    ii1(iii1)=iif; jj1(iii1)=iipzm; vv1(iii1)=1i/2;iii1=iii1+1;
    
%     %1/R parts
    g(iif)=g(iif)-1/R*(dfn-dfs)/(zn-zs);
    ii(iii)=iif; jj(iii)=iif; vv(iii)=1/R*1/(zpp-zp)/(zn-zs); iii=iii+1;
    ii(iii)=iif; jj(iii)=iif; vv(iii)=1/R*1/(zp-zpm)/(zn-zs);iii=iii+1;
    ii(iii)=iif; jj(iii)=iifzp; vv(iii)=1/R*-1/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iif; jj(iii)=iifzm; vv(iii)=1/R*-1/(zp-zpm)/(zn-zs);iii=iii+1;
    
%     g(iif)=g(iif)+1/R*u(iif)*U(iif)+1/R*U(iih)*(fn-fs)/(zn-zs)-2/R*U(iig)*u(iig)-2/R*u(iig)+u(iif)*bbar^2/R;
    g(iif)=g(iif)+1/R*u(iif)*U(iif)+1/R*U(iih)*(fn-fs)/(zn-zs)-2/R*U(iig)*u(iig)+u(iif)*bbar^2/R;
    
    ii(iii)=iif; jj(iii)=iif; vv(iii)=1/R*U(iif);iii=iii+1;
    
    ii(iii)=iif; jj(iii)=iif; vv(iii)=1/R*U(iih)*((zpp-zn))/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iif; jj(iii)=iifzp; vv(iii)=1/R*U(iih)*((zn-zp))/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iif; jj(iii)=iifzm; vv(iii)=-1/R*U(iih)*((zp-zs))/(zp-zpm)/(zn-zs);iii=iii+1;
    ii(iii)=iif; jj(iii)=iif; vv(iii)=-1/R*U(iih)*((zs-zpm))/(zp-zpm)/(zn-zs);iii=iii+1;
    
    ii(iii)=iif; jj(iii)=iig; vv(iii)=-2/R*U(iig);iii=iii+1;
%     ii(iii)=iif; jj(iii)=iig; vv(iii)=-2/R;iii=iii+1;
    ii(iii)=iif; jj(iii)=iif; vv(iii)=bbar^2/R;iii=iii+1;
    
    
    g(iif)=g(iif)+u(iif)*alpha^2/R;
    ii2(iii2)=iif; jj2(iii2)=iif; vv2(iii2)=1/R;iii2=iii2+1;
% 
%     
%     
%     %y mom
    gn=(u(iig)*(zpp-zn)+u(iigzp)*(zn-zp))/(zpp-zp);
    gs=(u(iigzm)*(zp-zs)+u(iig)*(zs-zpm))/(zp-zpm);
    Gn=(U(iig)*(zpp-zn)+U(iigzp)*(zn-zp))/(zpp-zp);
    Gs=(U(iigzm)*(zp-zs)+U(iig)*(zs-zpm))/(zp-zpm);
    pp=(u(iip)+u(iipzm))/2;
%     
        dgn=(u(iigzp)-u(iig))/(zpp-zp);
        dgs=(u(iig)-u(iigzm))/(zp-zpm);
    
    g(iig)=u(iig)*(1i*bbar*U(iig)-1i*omega)+u(iih)*(Gn-Gs)/(zn-zs)+1i*bbar*pp;
    g(iig)=g(iig)+u(iig)*1i*alpha*U(iif);
%     
%     %
    iio(iiio)=iig; jjo(iiio)=iig; vvo(iiio)=-1i; iiio=iiio+1;

%     ii(iii)=iig; jj(iii)=iig; vv(iii)=(1i*bbar*U(iig)-1i*omega);iii=iii+1;
    ii(iii)=iig; jj(iii)=iig; vv(iii)=(1i*bbar*U(iig));iii=iii+1;
    ii(iii)=iig; jj(iii)=iih; vv(iii)=(Gn-Gs)/(zn-zs);iii=iii+1;
    ii(iii)=iig; jj(iii)=iip; vv(iii)=bbar*1i/2;iii=iii+1;
    ii(iii)=iig; jj(iii)=iipzm; vv(iii)=bbar*1i/2;iii=iii+1;
    
    ii1(iii1)=iig; jj1(iii1)=iig; vv1(iii1)=1i*U(iif);iii1=iii1+1;
    
%     %1/R terms
    g(iig)=g(iig)-1/R*(dgn-dgs)/(zn-zs);
    ii(iii)=iig; jj(iii)=iig; vv(iii)=1/R*1/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iig; jj(iii)=iig; vv(iii)=1/R*1/(zp-zpm)/(zn-zs);iii=iii+1;
    ii(iii)=iig; jj(iii)=iigzp; vv(iii)=1/R*-1/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iig; jj(iii)=iigzm; vv(iii)=1/R*-1/(zp-zpm)/(zn-zs);iii=iii+1;
%     
%     g(iig)=g(iig)+2/R*u(iif)*U(iig)+1/R*U(iif)*u(iig)+1/R*U(iih)*(gn-gs)/(zn-zs)+2/R*u(iif)+1/R*bbar^2*u(iig);
    g(iig)=g(iig)+2/R*u(iif)*U(iig)+1/R*U(iif)*u(iig)+1/R*U(iih)*(gn-gs)/(zn-zs)+1/R*bbar^2*u(iig);
    
    ii(iii)=iig; jj(iii)=iif; vv(iii)=2/R*U(iig);iii=iii+1;
    ii(iii)=iig; jj(iii)=iig; vv(iii)=1/R*U(iif);iii=iii+1;
    
    ii(iii)=iig; jj(iii)=iig; vv(iii)=1/R*U(iih)*((zpp-zn))/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iig; jj(iii)=iigzp; vv(iii)=1/R*U(iih)*((zn-zp))/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iig; jj(iii)=iig; vv(iii)=-1/R*U(iih)*((zs-zpm))/(zp-zpm)/(zn-zs);iii=iii+1;
    ii(iii)=iig; jj(iii)=iigzm; vv(iii)=-1/R*U(iih)*((zp-zs))/(zp-zpm)/(zn-zs);iii=iii+1;
    
%     ii(iii)=iig; jj(iii)=iif; vv(iii)=2/R;iii=iii+1;
    ii(iii)=iig; jj(iii)=iig; vv(iii)=1/R*bbar^2;iii=iii+1;
    
    

    g(iig)=g(iig)+1/R*alpha^2*u(iig);
    ii2(iii2)=iig; jj2(iii2)=iig; vv2(iii2)=1/R;iii2=iii2+1;
%     
%     
%     
%     %zmom
    Hn=(U(iih)*(zpp-zn)+U(iihzp)*(zn-zp))/(zpp-zp);
    Hs=(U(iihzm)*(zp-zs)+U(iih)*(zs-zpm))/(zp-zpm);
    
    hn=(u(iih)*(zpp-zn)+u(iihzp)*(zn-zp))/(zpp-zp);
    hs=(u(iihzm)*(zp-zs)+u(iih)*(zs-zpm))/(zp-zpm);
    
    dhn=(u(iihzp)-u(iih))/(zpp-zp);
    dhs=(u(iih)-u(iihzm))/(zp-zpm);
    
    g(iih)=u(iih)*(1i*bbar*U(iig)-1i*omega)+(u(iip)-u(iipzm))/(zn-zs);
    
        iio(iiio)=iih; jjo(iiio)=iih; vvo(iiio)=-1i; iiio=iiio+1;

    
%     ii(iii)=iih; jj(iii)=iih; vv(iii)=(1i*bbar*U(iig)-1i*omega); iii=iii+1;
    ii(iii)=iih; jj(iii)=iih; vv(iii)=(1i*bbar*U(iig)); iii=iii+1;
    ii(iii)=iih; jj(iii)=iip; vv(iii)=1/(zn-zs);iii=iii+1;
    ii(iii)=iih; jj(iii)=iipzm; vv(iii)=-1/(zn-zs);iii=iii+1;
    
    g(iih)=g(iih)+u(iih)*1i*alpha*U(iif);
    ii1(iii1)=iih; jj1(iii1)=iih; vv1(iii1)=1i*U(iif);iii1=iii1+1;
    
%     %1/R parts
    g(iih)=g(iih)+u(iih)*alpha^2/R;
    ii2(iii2)=iih; jj2(iii2)=iih; vv2(iii2)= 1/R;iii2=iii2+1;
%     
    g(iih)=g(iih)+u(iih)*(bbar.^2/R+1/R*(Hn-Hs)/(zn-zs));
    ii(iii)=iih; jj(iii)=iih; vv(iii)=(bbar.^2/R+1/R*(Hn-Hs)/(zn-zs));iii=iii+1;
%     
    g(iih)=g(iih)-1/R*(dhn-dhs)/(zn-zs)+1/R*U(iih)*(hn-hs)/(zn-zs);
    
    ii(iii)=iih; jj(iii)=iihzp; vv(iii)=-1/R*1/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iih; jj(iii)=iih; vv(iii)=-1/R*-1/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iih; jj(iii)=iih; vv(iii)=-1*-1/R*1/(zp-zpm)/(zn-zs);iii=iii+1;
    ii(iii)=iih; jj(iii)=iihzm; vv(iii)=-1*-1/R*-1/(zp-zpm)/(zn-zs);iii=iii+1;

    ii(iii)=iih; jj(iii)=iihzp; vv(iii)=1/R*U(iih)*((zn-zp))/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iih; jj(iii)=iih; vv(iii)=1/R*U(iih)*((zpp-zn))/(zpp-zp)/(zn-zs);iii=iii+1;
    ii(iii)=iih; jj(iii)=iihzm; vv(iii)=-1/R*U(iih)*((zp-zs))/(zp-zpm)/(zn-zs);iii=iii+1;
    ii(iii)=iih; jj(iii)=iih; vv(iii)=-1/R*U(iih)*((zs-zpm))/(zp-zpm)/(zn-zs);iii=iii+1;
%     
end
% 
% %bc

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
% % 
% % 
% % i=np; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
% % g(iif)=u(iif); g(iig)=u(iig); 
% % 
% % % g(iih)=u(5); ii(iii)=iih]; jj(iii)=5]; vv(iii)=1];
% % g(iih)=u(iih); ii(iii)=iih]; jj(iii)=iih]; vv(iii)=1]; %artificial contion but helped a bit, maybe we should fix p to 0 at 0,inf also
% % ii(iii)=iif]; jj(iii)=iif]; vv(iii)=1];
% % ii(iii)=iig]; jj(iii)=iig]; vv(iii)=1];
% 
i=length(zw); iip=(i-1)*4+1; iif=(i-1)*4+2; iig=(i-1)*4+3; iih=(i-1)*4+4;
iifzm=iif-4; iigzm=iig-4; iihzm=iih-4;
g(iif)=(u(iif)+u(iifzm))/2; g(iig)=(u(iig)+u(iigzm))/2;  g(iip)=u(iip)-u(iip-4);

ii(iii)=iif; jj(iii)=iif; vv(iii)=1/2;iii=iii+1;
ii(iii)=iif; jj(iii)=iifzm; vv(iii)=1/2;iii=iii+1;
ii(iii)=iig; jj(iii)=iig; vv(iii)=1/2;iii=iii+1;
ii(iii)=iig; jj(iii)=iigzm; vv(iii)=1/2;iii=iii+1;
ii(iii)=iip; jj(iii)=iip; vv(iii)=1;iii=iii+1;
ii(iii)=iip; jj(iii)=iip-4; vv(iii)=-1;iii=iii+1;
% 
% 
% g(iih)=u(1); %pressure pert at wall = 0 %this generates not full rank jac
% ii(iii)=iih; jj(iii)=1; vv(iii)=1;iii=iii+1;

g(iih)=(u(iih)+u(iihzm))/2; %h pert at inf somehow artificial
ii(iii)=iih; jj(iii)=iih; vv(iii)=1/2;iii=iii+1;
ii(iii)=iih; jj(iii)=iihzm; vv(iii)=1/2;iii=iii+1;


% i=1; iip=(i-1)*4+1; g(iip)=u(iip); ii(iii)=iip; jj(iii)=iip; vv(iii)=1; iii=iii+1;
% i=length(zw)-1; iip=(i-1)*4+1; g(iip)=u(iip); ii(iii)=iip; jj(iii)=iip; vv(iii)=1; iii=iii+1;

jac0=sparse(ii,jj,vv);
jac1=sparse(ii1,jj1,vv1);
jac2=sparse(ii2,jj2,vv2);
jacom=sparse(iio,jjo,vvo);
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