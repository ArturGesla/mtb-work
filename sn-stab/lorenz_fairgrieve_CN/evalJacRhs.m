
for ip=1:np-1
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

    T=u(neq*np+1); dt=T/(np-1); ds=1/(np-1);

    x=(u(ixp)+u(ix))/2;
    y=(u(iyp)+u(iy))/2;
    z=(u(izp)+u(iz))/2;
    
    g(ix)=T*(sigma*(y-x))           +(u(ix)-u(ixp))/ds;
    g(iy)=T*(r*x-y-x*z)     +(u(iy)-u(iyp))/ds;
    g(iz)=T*(x*y-b*z)           +(u(iz)-u(izp))/ds;

    
    J(ix,ix)=T*(sigma*(-1/2))+1/ds;         J(ix,ixp)=T*(sigma*(-1/2))+(-1)/ds;
    J(ix,iy)=T*(sigma*(1/2));               J(ix,iyp)=T*(sigma*(1/2));
    J(ix,neq*np+1)=(sigma*(y-x));

    J(iy,iy)=T*(-1/2)+1/ds;                 J(iy,iyp)=T*(-1/2)-1/ds;
    J(iy,ix)=T*(r/2-z/2);                   J(iy,ixp)=T*(r/2-z/2);
    J(iy,iz)=T*(-x/2);                      J(iy,izp)=T*(-x/2);
    J(iy,neq*np+1)=(r*x-y-x*z);

    J(iz,iz)=T*(-b)/2+1/ds;                   J(iz,izp)=T*(-b)/2-1/ds;
    J(iz,ix)=T*(y/2);                           J(iz,ixp)=T*(y/2);
    J(iz,iy)=T*(x/2);                       J(iz,iyp)=T*(x/2);
    J(iz,neq*np+1)=(x*y-b*z);
    
    
end
    
ip=np;
ix=ip*neq-2;
iy=ip*neq-1;
iz=ip*neq;
        
J(ix,ix)=1; J(ix,1)=-1; g(ix)=u(ix)-u(1);
J(iy,iy)=1; J(iy,2)=-1; g(iy)=u(iy)-u(2);
J(iz,iz)=1; J(iz,3)=-1; g(iz)=u(iz)-u(3);

%
    % last row
%     g(iz+1)=(u(ix)-u(ixpp))/ds/2-60.9946;
%     J(iz+1,ix)=1/ds/2;
%     J(iz+1,ixpp)=-1/ds/2;
    
%     g(iz+1)=(u(ixm)-u(ixp))/ds/2-60.9946;
% 

%% phase cond
% ip=phaseIndex;
% % if (ip==1 && ip==np) warning("ip lols - fatal"); end
% ix=ip*neq-2;
%         ixp=mod(ip*neq-2+neq-1,neq*np)+1;
%         ixm=mod(ip*neq-2-neq-1,neq*np)+1-3*(ip==1);
%        
% 
% g(iz+1)=(u(ixm)-u(ixp))/ds/2-derX;
% % g(iz+1)=(u(ixm)-u(ixp))/ds/2+13.212372967119705;
%     J(iz+1,ixm)=1/ds/2;
%     J(iz+1,ixp)=-1/ds/2;

%% new pc
ip=1;
ix=ip*neq-2; ixp=ix+neq; ixm=neq*np-2-3;
iy=ip*neq-1; iyp=iy+neq; iym=neq*np-1-3;
iz=ip*neq; izp=iz+neq; izm=neq*np-3;


% g(iz+1)=(x0p-x0m)/ds/2*(u(ix)-x0);
J(np*neq+1,ix)=(u(ixp)-u(ixm))/ds/2;
J(np*neq+1,iy)=(u(iyp)-u(iym))/ds/2;
J(np*neq+1,iz)=(u(izp)-u(izm))/ds/2;