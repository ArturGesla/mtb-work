
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

%     [a*y(1)+b*y(2)+y(1)*y(3)+c2*y(1)*(1-y(3)^2);
%                     c*y(1)+d*y(2)+y(2)*y(3)+c2*y(2)*(1-y(3)^2);
%                     e*y(3)-(y(1)*y(1)+y(2)*y(2)+y(3)*y(3))];% Anonymous Function
    
    g(ix)=T*(a*x+b*y+x*z+c2*x*(1-z^2))                +(u(ix)-u(ixp))/ds;
    g(iy)=T*(c*x+d*y+y*z+c2*y*(1-z^2))                +(u(iy)-u(iyp))/ds;
    g(iz)=T*(e*z-(x*x+y*y+z*z))                        +(u(iz)-u(izp))/ds;

    
    J(ix,ix)=T*(a*1/2+1/2*z+c2*1/2*(1-z^2))+(1)/ds;         J(ix,ixp)=T*(a*1/2+1/2*z+c2*1/2*(1-z^2))+(-1)/ds;
    J(ix,iy)=T*(b*1/2);                                     J(ix,iyp)=T*(b*1/2);
    J(ix,iz)=T*(x*1/2+c2*x*(-2*z*1/2));                     J(ix,izp)=T*(x*1/2+c2*x*(-2*z*1/2));
    J(ix,neq*np+1)=(a*x+b*y+x*z+c2*x*(1-z^2));
% 
    J(iy,iy)=T*(d*1/2+z*1/2+c2*1/2*(1-z^2))+1/ds;           J(iy,iyp)=T*(d*1/2+z*1/2+c2*1/2*(1-z^2))-1/ds;
    J(iy,ix)=T*(1/2*c);                                     J(iy,ixp)=T*(1/2*c);
    J(iy,iz)=T*(1/2*y+c2*y*(-2*z*1/2));                     J(iy,izp)=T*(1/2*y+c2*y*(-2*z*1/2));
    J(iy,neq*np+1)=(c*x+d*y+y*z+c2*y*(1-z^2));

    J(iz,iz)=T*(e*1/2-(2*z*1/2))+(1)/ds;                    J(iz,izp)=T*(e*1/2-(2*z*1/2))+(-1)/ds;
    J(iz,ix)=T*(-(2*x*1/2));                                J(iz,ixp)=T*(-(2*x*1/2));
    J(iz,iy)=T*(-(2*y*1/2));                                       J(iz,iyp)=T*(-(2*y*1/2));
    J(iz,neq*np+1)=(e*z-(x*x+y*y+z*z));
    
    
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
ip=phaseIndex;
if (ip==1 && ip==np) warning("ip lols - fatal"); end
ix=ip*neq-2;
        ixp=mod(ip*neq-2+neq-1,neq*np)+1;
        ixm=mod(ip*neq-2-neq-1,neq*np)+1;
       

g(iz+1)=(u(ixm)-u(ixp))/ds/2-derX;
% g(iz+1)=(u(ixm)-u(ixp))/ds/2+13.212372967119705;
    J(iz+1,ixm)=1/ds/2;
    J(iz+1,ixp)=-1/ds/2;
