function [rhs,jac,B]=calculateJacAndRhs(zc,zw,u,Re,ra,kr,L)
%%

jac=sparse([]);
B=sparse([]);
rhs=[];

Nz=length(zc);
%%


for iz=2:Nz-1
    indices;
%     
%     %conti
%     mult=1;
%     rhs(iip)=2*u(iif)+mult*(u(iih)-u(iihzm))/hz;
%     
%     jac(iip,iif)=2+ra*1i*kr;
%     jac(iip,iih)=mult*1/hz;
%     jac(iip,iihzm)=mult*-1/hz;
    
%     rhs(iip)=2*u(iif)+(u(iih)-u(iihzm))/hz;
%     
%     jac(iip,iif)=2;
%     jac(iip,iih)=1/hz;
%     jac(iip,iihzm)=-1/hz;
    
    %r mom
    rhs(iif)=u(iif)^2-u(iig)^2+u(iih)*(fs-fn)/hz-1/Re*(azup*u(iifzp)+azum*u(iifzm)-u(iif)*(azum+azup))+u(iik);
    
    jac(iif,iif)=2*u(iif)-u(iih)*((zw(iz-1)-zc(iz-1))*1)/(zc(iz)-zc(iz-1))/hz+u(iih)*((zc(iz+1)-zw(iz))*1)/(zc(iz+1)-zc(iz))/hz-1/Re*(-1*(azum+azup))+ra*u(iif)*1i*kr-L^2/Re*(3*1i*kr/ra-kr^2);
    jac(iif,iifzm)=-u(iih)*((zc(iz)-zw(iz-1))*1)/(zc(iz)-zc(iz-1))/hz    -1/Re*(azum);
    jac(iif,iifzp)=u(iih)*((zw(iz)-zc(iz))*1)/(zc(iz+1)-zc(iz))/hz      -1/Re*(azup);
    
    jac(iif,iig)=-2*u(iig);
    jac(iif,iih)=(fs-fn)/hz;
    
%     if kr==0
if kr==0 || kr~=1
    jac(iif,iik)=1;%+1/2*ra*1i*kr;
    end
    jac(iif,iip)=L^2/ra*1i*kr;
 
    B(iif,iif)=-1;
    %theta mom
    rhs(iig)=2*u(iif)*u(iig)+u(iih)*(gs-gn)/hz-1/Re*(azvp*u(iigzp)+azvm*u(iigzm)-u(iig)*(azvm+azvp));
    
    jac(iig,iig)=2*u(iif)-u(iih)*((zw(iz-1)-zc(iz-1))*1)/(zc(iz)-zc(iz-1))/hz+u(iih)*((zc(iz+1)-zw(iz))*1)/(zc(iz+1)-zc(iz))/hz-1/Re*(-1*(azvm+azvp))+ra*u(iif)*1i*kr-L^2/Re*(3*1i*kr/ra-kr^2);
    jac(iig,iigzm)=-u(iih)*((zc(iz)-zw(iz-1))*1)/(zc(iz)-zc(iz-1))/hz    -1/Re*(azvm);
    jac(iig,iigzp)=u(iih)*((zw(iz)-zc(iz))*1)/(zc(iz+1)-zc(iz))/hz      -1/Re*(azvp);
    
    jac(iig,iif)=2*u(iig);
    jac(iig,iih)=(gs-gn)/hz;
    
    B(iig,iig)=-1;
end

for iz=2:Nz-2
    indices;
    
        
    %conti
    mult=1;
    rhs(iip)=2*u(iif)+mult*(u(iih)-u(iihzm))/hz;
    
    jac(iip,iif)=2+ra*1i*kr;
    jac(iip,iih)=mult*1/hz;
    jac(iip,iihzm)=mult*-1/hz;


    %z mom
    rhs(iih)=u(iih)*(hs-hn)/hzw+(ps-pn)/hzw-1/Re*(azwp*u(iihzp)+azwm*u(iihzm)-u(iih)*(azwm+azwp));
    
    jac(iih,iih)=(hs-hn)/hzw+u(iih)*((zw(iz+1)-zc(iz+1))*1)/(zw(iz+1)-zw(iz))/hzw-u(iih)*((zc(iz)-zw(iz-1))*1)/(zw(iz)-zw(iz-1))/hzw-1/Re*(-1*(azwm+azwp))+ra*u(iif)*1i*kr-L^2/Re*(1i*kr/ra-kr^2);
    jac(iih,iihzm)=-u(iih)*((zw(iz)-zc(iz))*1)/(zw(iz)-zw(iz-1))/hzw    -1/Re*(azwm);
    jac(iih,iihzp)=u(iih)*((zc(iz+1)-zw(iz))*1)/(zw(iz+1)-zw(iz))/hzw      -1/Re*(azwp);
    
    jac(iih,iip)=-1/hzw;
    jac(iih,iipzp)=1/hzw; 
    
    B(iih,iih)=-1;
end

%last conti
%conti
iz=Nz-1;
indices;
rhs(iip)=u(iip);
jac(iip,iip)=1;

%boundaries
%top
iz=1;
indices;
rhs(iip)=u(iip);
rhs(iif)=u(iif)+u(iifzp);
rhs(iig)=u(iig)+u(iigzp);
rhs(iih)=u(iih);

jac(iip,iip)=1;
jac(iif,iif)=1;
jac(iif,iifzp)=1;
jac(iig,iig)=1;
jac(iig,iigzp)=1;
jac(iih,iih)=1;

%bottom
iz=Nz;
indices;
rhs(iip)=u(iip);
rhs(iif)=u(iif)+u(iifzm);
rhs(iig)=u(iig)+u(iigzm)-1*2;
rhs(iih)=u(iih);
rhs(iihzm)=u(iihzm);

jac(iip,iip)=1;
jac(iif,iif)=1;
jac(iif,iifzm)=1;
jac(iig,iig)=1;
jac(iig,iigzm)=1;
jac(iih,iih)=1;
jac(iihzm,iihzm)=1;
B(iih,iih)=0;

if kr==0 || kr~=1
%last
iz=Nz-1;
indices;


%conti
    mult=1;
    rhs(iik)=2*u(iif)+mult*(u(iih)-u(iihzm))/hz;
    
    jac(iik,iif)=2+ra*1i*kr;
    jac(iik,iih)=mult*1/hz;
    jac(iik,iihzm)=mult*-1/hz;
    
%     
% indices;
% rhs(iik)=u(iipzp);
% jac(iik,iipzp)=1;

B(iik,iik)=0;
end
end

