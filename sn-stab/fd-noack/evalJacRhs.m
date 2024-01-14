% f=[mu*u-v*g-w*u;
%     mu*v+u/g-v*w;
%     -w+u^2+v^2*g*g];

for ip=1:np-1
        ix=ip*neq-2;
        iy=ip*neq-1;
        iz=ip*neq;

        ixp=mod(ip*neq-2+neq-1,neq*np)+1;
%         ixpp=mod(ip*neq-2+neq+neq-1,neq*np)+1;
        iyp=mod(ip*neq-1+neq-1,neq*np)+1;
        izp=mod(ip*neq-0+neq-1,neq*np)+1;
%         ixm=mod(ip*neq-2-neq-1,neq*np)+1-neq*((ip*neq-2-neq-1)<0);
%         iym=mod(ip*neq-1-neq-1,neq*np)+1-neq*((ip*neq-1-neq-1)<0);
%         izm=mod(ip*neq-0-neq-1,neq*np)+1-neq*((ip*neq-0-neq-1)<0);

    T=u(neq*np+1); 
%     dt=T/(np-1);
    ds=1/(np-1);

    x=(u(ixp)+u(ix))/2;
    y=(u(iyp)+u(iy))/2;
    z=(u(izp)+u(iz))/2;
    
    g(ix)=T*(mu*x-y*gm-z*x)           +(u(ix)-u(ixp))/ds;
    g(iy)=T*(mu*y+x/gm-y*z)     +(u(iy)-u(iyp))/ds;
    g(iz)=T*(-z+x*x+y*y*gm*gm)           +(u(iz)-u(izp))/ds;

    
    J(ix,ix)=T*(mu*(1/2)-z/2)+1/ds;         J(ix,ixp)=T*(mu*(1/2)-z/2)+(-1)/ds;
    J(ix,iy)=T*(-gm/2);               J(ix,iyp)=T*(-gm/2);
    J(ix,iz)=T*(-x/2);               J(ix,izp)=T*(-x/2);
    J(ix,neq*np+1)=(mu*x-y*gm-z*x);
% 
    J(iy,iy)=T*(mu-z)/2+1/ds;                 J(iy,iyp)=T*(mu-z)/2-1/ds;
    J(iy,ix)=T*(1/gm/2);                   J(iy,ixp)=T*(1/gm/2);
    J(iy,iz)=T*(-y/2);                      J(iy,izp)=T*(-y/2);
    J(iy,neq*np+1)=(mu*y+x/gm-y*z);
% 
    J(iz,iz)=T*(-1)/2+1/ds;                   J(iz,izp)=T*(-1)/2-1/ds;
    J(iz,ix)=T*(x);                           J(iz,ixp)=T*(x);
    J(iz,iy)=T*(y*gm*gm);                       J(iz,iyp)=T*(y*gm*gm);
    J(iz,neq*np+1)=(-z+x*x+y*y*gm*gm);
%     
    
end
%     
ip=np;
ix=ip*neq-2;
iy=ip*neq-1;
iz=ip*neq;
        
J(ix,ix)=1; J(ix,1)=-1; g(ix)=u(ix)-u(1);
J(iy,iy)=1; J(iy,2)=-1; g(iy)=u(iy)-u(2);
J(iz,iz)=1; J(iz,3)=-1; g(iz)=u(iz)-u(3);
% 
% %
%     % last row
% %     g(iz+1)=(u(ix)-u(ixpp))/ds/2-60.9946;
% %     J(iz+1,ix)=1/ds/2;
% %     J(iz+1,ixpp)=-1/ds/2;
%     
% %     g(iz+1)=(u(ixm)-u(ixp))/ds/2-60.9946;
% % 
% ip=phaseIndex;
% if (ip==1 && ip==np) warning("ip lols - fatal"); end
ip=1;
ix=ip*neq-2;
        ixp=mod(ip*neq-2+neq-1,neq*np)+1;
%         ixm=mod(ip*neq-2-neq-1,neq*np)+1;
%        
% 
g(iz+1)=(u(ix)-u(ixp))/ds;
% % g(iz+1)=(u(ixm)-u(ixp))/ds/2+13.212372967119705;
    J(iz+1,ix)=1/ds;
    J(iz+1,ixp)=-1/ds;
