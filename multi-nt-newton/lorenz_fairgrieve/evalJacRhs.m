
for ip=2:np-1
        ix=ip*neq-2;
        iy=ip*neq-1;
        iz=ip*neq;

        ixp=mod(ip*neq-2+neq-1,neq*np)+1;
%         ixpp=mod(ip*neq-2+neq+neq-1,neq*np)+1;
        iyp=mod(ip*neq-1+neq-1,neq*np)+1;
        izp=mod(ip*neq-0+neq-1,neq*np)+1;
        ixm=mod(ip*neq-2-neq-1,neq*np)+1-neq*((ip*neq-2-neq-1)<0);
        iym=mod(ip*neq-1-neq-1,neq*np)+1-neq*((ip*neq-1-neq-1)<0);
        izm=mod(ip*neq-0-neq-1,neq*np)+1-neq*((ip*neq-0-neq-1)<0);

%     x=u(ix); y=u(iy); z=u(iz);% r=sqrt(x^2+y^2);
% 
%     xNext=u(ixp); yNext=u(iyp); zNext=u(izp);
%     xPrev=u(ixm); yPrev=u(iym); zPrev=u(izm);

%     xNextNext=u(mod(ip*2-1+neq+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);

    T=u(neq*np+1); dt=T/(np-2); ds=1/(np-2);

    g(ix)=T*(sigma*(u(iy)-u(ix)))           +(u(ixm)-u(ixp))/ds/2;
    g(iy)=T*(r*u(ix)-u(iy)-u(ix)*u(iz))     +(u(iym)-u(iyp))/ds/2;
    g(iz)=T*(u(ix)*u(iy)-b*u(iz))           +(u(izm)-u(izp))/ds/2;

    gstab(ix)=T*(sigma*(u(iy)-u(ix)))      ;
    gstab(iy)=T*(r*u(ix)-u(iy)-u(ix)*u(iz))   ;
    gstab(iz)=T*(u(ix)*u(iy)-b*u(iz))           ;



    
    
    J(ix,ix)=T*(sigma*(-1));
    J(ix,iy)=T*(sigma*(1));
    J(ix,ixp)=(-1)/ds/2;
    J(ix,ixm)=(1)/ds/2;
    J(ix,neq*np+1)=(sigma*(u(iy)-u(ix)));

    J(iy,iy)=T*(-1);
    J(iy,ix)=T*(r-u(iz));
    J(iy,iz)=T*(-u(ix));
    J(iy,iyp)=(-1)/ds/2;
    J(iy,iym)=(1)/ds/2;
    J(iy,neq*np+1)=(r*u(ix)-u(iy)-u(ix)*u(iz));

    J(iz,iz)=T*(-b);
    J(iz,ix)=T*(u(iy));
    J(iz,iy)=T*(u(ix));
    J(iz,izp)=(-1)/ds/2;
    J(iz,izm)=(1)/ds/2;
    J(iz,neq*np+1)=(u(ix)*u(iy)-b*u(iz));
    
    %jstab
    
     Jstab(ix,ix)=T*(sigma*(-1));
    Jstab(ix,iy)=T*(sigma*(1));
%     Jstab(ix,ixp)=(-1)/ds/2;
%     Jstab(ix,ixm)=(1)/ds/2;
%     J(ix,neq*np+1)=(sigma*(u(iy)-u(ix)));

    Jstab(iy,iy)=T*(-1);
    Jstab(iy,ix)=T*(r-u(iz));
    Jstab(iy,iz)=T*(-u(ix));
%     Jstab(iy,iyp)=(-1)/ds/2;
%     Jstab(iy,iym)=(1)/ds/2;
%     J(iy,neq*np+1)=(r*u(ix)-u(iy)-u(ix)*u(iz));

    Jstab(iz,iz)=T*(-b);
    Jstab(iz,ix)=T*(u(iy));
    Jstab(iz,iy)=T*(u(ix));
%     Jstab(iz,izp)=(-1)/ds/2;
%     Jstab(iz,izm)=(1)/ds/2;
%     J(iz,neq*np+1)=(u(ix)*u(iy)-b*u(iz));
  
    end
%% almost last row - effectively the repeated point
ip=np-1;

ix=ip*neq-2;
iy=ip*neq-1;
iz=ip*neq;
        
J(1,ix)=1; J(1,1)=-1; g(1)=u(ix)-u(1);
J(2,iy)=1; J(2,2)=-1; g(2)=u(iy)-u(2);
J(3,iz)=1; J(3,3)=-1; g(3)=u(iz)-u(3);

ip=np;

ix=ip*neq-2;
iy=ip*neq-1;
iz=ip*neq;
        
J(ix,ix)=1; J(ix,1+neq)=-1; g(ix)=u(ix)-u(1+neq);
J(iy,iy)=1; J(iy,2+neq)=-1; g(iy)=u(iy)-u(2+neq);
J(iz,iz)=1; J(iz,3+neq)=-1; g(iz)=u(iz)-u(3+neq);

%%
    % last row
%     g(iz+1)=(u(ix)-u(ixpp))/ds/2-60.9946;
%     J(iz+1,ix)=1/ds/2;
%     J(iz+1,ixpp)=-1/ds/2;
    
%     g(iz+1)=(u(ixm)-u(ixp))/ds/2-60.9946;
% 
ip=phaseIndex;
if (ip==1 && ip==np) warning("ip lols - fatal"); end
ix=ip*neq-2;
        ixp=mod(ip*neq-2+neq-1,neq*np)+1;
        ixm=mod(ip*neq-2-neq-1,neq*np)+1;
       

g(iz+1)=(u(ixm)-u(ixp))/ds/2-derX;
% g(iz+1)=(u(ixm)-u(ixp))/ds/2+13.212372967119705;
    J(iz+1,ixm)=1/ds/2;
    J(iz+1,ixp)=-1/ds/2;
