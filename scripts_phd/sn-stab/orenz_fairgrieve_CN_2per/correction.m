 for i=1:8
    for ip=1:np
        ix=ip*neq-2;
        iy=ip*neq-1;
        iz=ip*neq;

        ixp=mod(ip*neq-2+neq-1,neq*np)+1;
        ixpp=mod(ip*neq-2+neq+neq-1,neq*np)+1;
        iyp=mod(ip*neq-1+neq-1,neq*np)+1;
        izp=mod(ip*neq-0+neq-1,neq*np)+1;
        ixm=mod(ip*neq-2-neq-1,neq*np)+1;
        iym=mod(ip*neq-1-neq-1,neq*np)+1;
        izm=mod(ip*neq-0-neq-1,neq*np)+1;

%     x=u(ix); y=u(iy); z=u(iz);% r=sqrt(x^2+y^2);
% 
%     xNext=u(ixp); yNext=u(iyp); zNext=u(izp);
%     xPrev=u(ixm); yPrev=u(iym); zPrev=u(izm);

%     xNextNext=u(mod(ip*2-1+neq+neq-1,neq*np)+1); yNext=u(mod(ip*2+neq-1,neq*np)+1);

    T=u(neq*np+1); dt=T/np; ds=1/np;

    g(ix)=T*(sigma*(u(iy)-u(ix)))           +(u(ixm)-u(ixp))/ds/2;
    g(iy)=T*(r*u(ix)-u(iy)-u(ix)*u(iz))     +(u(iym)-u(iyp))/ds/2;
    g(iz)=T*(u(ix)*u(iy)-b*u(iz))           +(u(izm)-u(izp))/ds/2;

    
    
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
    
     %Jlam
    Jlam(iy)=T*u(ix);
  
    end

    % last row
    g(iz+1)=(u(ix)-u(ixpp))/ds/2-60.9946;
    J(iz+1,ix)=1/ds/2;
    J(iz+1,ixpp)=-1/ds/2;
    
   

    %
    du=-sparse(J)\g;
   
    fprintf('%s\n',"iter: "+num2str(i)+" res: "+num2str(norm(du)))
    u=u+du;
    uM=[uM,u];
 end