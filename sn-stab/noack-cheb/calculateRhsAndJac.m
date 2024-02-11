function [g,jac]=calculateRhsAndJac(neq,nt,u,r,om1)
%%

s=10;
b=8/3;
% r=24;
mu=r.^2;
ii=[]; jj=[]; vv=[];

g=u*0;

%u=[xr0,yr0,zr0,xr1,...,xi0,yi0,zi0,...]

% Noack system
% f=[mu*u-v*g-w*u;
%     mu*v+u/g-v*w;
%     -w+u^2+v^2*g*g];

lastid=nt*neq+1;
om=u(lastid);
% for ik=0:nt-1%2*nt-1
    for ik=0:nt-1-1% last one is evp
    ix=ik*neq+1; iy=ix+1; iz=iy+1;
    
%     nlxz=0;
%     nlxy=0;
    
        %non linear
    for ikl=ik:-1:0
        ikr=ik-ikl;
        
        ixr=abs(ikr)*neq+1;  ixl=abs(ikl)*neq+1; iyr=abs(ikr)*neq+2;  iyl=abs(ikl)*neq+2; izr=abs(ikr)*neq+3;  izl=abs(ikl)*neq+3; 
        
        g(ix)=g(ix)-u(ixl)*u(izr)/2;
        g(iy)=g(iy)-u(izl)*u(iyr)/2;
        g(iz)=g(iz)+u(ixl)*u(ixr)/2+u(iyl)*u(iyr)/2;

        ii(end+1)=ix; jj(end+1)=ixl; vv(end+1)=-u(izr)/2;
         ii(end+1)=ix; jj(end+1)=izr; vv(end+1)=-u(ixl)/2;
         
         ii(end+1)=iy; jj(end+1)=izl; vv(end+1)=-u(iyr)/2;
         ii(end+1)=iy; jj(end+1)=iyr; vv(end+1)=-u(izl)/2;
         
         ii(end+1)=iz; jj(end+1)=ixl; vv(end+1)=u(ixr)/2;
         ii(end+1)=iz; jj(end+1)=ixr; vv(end+1)=u(ixl)/2;
         ii(end+1)=iz; jj(end+1)=iyr; vv(end+1)=u(iyl)/2;
         ii(end+1)=iz; jj(end+1)=iyl; vv(end+1)=u(iyr)/2;
%          
    end
    for ikl=ik:nt-1
        ikr=ikl-ik;
        
        ixr=abs(ikr)*neq+1;  ixl=abs(ikl)*neq+1; iyr=abs(ikr)*neq+2;  iyl=abs(ikl)*neq+2; izr=abs(ikr)*neq+3;  izl=abs(ikl)*neq+3; 
        
        g(ix)=g(ix)-u(ixl)*u(izr)/2*(1+(ik~=0));
        g(iy)=g(iy)-u(izl)*u(iyr)/2*(1+(ik~=0));
        g(iz)=g(iz)+u(ixl)*u(ixr)/2*(1+(ik~=0))+u(iyl)*u(iyr)/2*(1+(ik~=0));

        ii(end+1)=ix; jj(end+1)=ixl; vv(end+1)=-u(izr)/2*(1+(ik~=0));
         ii(end+1)=ix; jj(end+1)=izr; vv(end+1)=-u(ixl)/2*(1+(ik~=0));
         
         ii(end+1)=iy; jj(end+1)=izl; vv(end+1)=-u(iyr)/2*(1+(ik~=0));
         ii(end+1)=iy; jj(end+1)=iyr; vv(end+1)=-u(izl)/2*(1+(ik~=0));
         
         ii(end+1)=iz; jj(end+1)=ixl; vv(end+1)=u(ixr)/2*(1+(ik~=0));
         ii(end+1)=iz; jj(end+1)=ixr; vv(end+1)=u(ixl)/2*(1+(ik~=0));
         ii(end+1)=iz; jj(end+1)=iyr; vv(end+1)=u(iyl)/2*(1+(ik~=0));
         ii(end+1)=iz; jj(end+1)=iyl; vv(end+1)=u(iyr)/2*(1+(ik~=0));
        
    end     
    
%linear terms
    g(ix)=g(ix)+mu*u(ix)-u(iy);
    g(iy)=g(iy)+mu*u(iy)+u(ix);
    g(iz)=g(iz)-u(iz);
    
    ii(end+1)=ix; jj(end+1)=iy; vv(end+1)=-1;
    ii(end+1)=ix; jj(end+1)=ix; vv(end+1)=mu;
    
    ii(end+1)=iy; jj(end+1)=ix; vv(end+1)=1;
    ii(end+1)=iy; jj(end+1)=iy; vv(end+1)=mu;
        
    ii(end+1)=iz; jj(end+1)=iz; vv(end+1)=-1;

%temporal term
    for ikl=ik+1:2:nt-1
    ixl=abs(ikl)*neq+1; iyl=abs(ikl)*neq+2; izl=abs(ikl)*neq+3; 
  
    g(ix)=g(ix)-u(ixl)*ikl*(2-(1*(ik==0)))*om/pi*1/om1;
    g(iy)=g(iy)-u(iyl)*ikl*(2-(1*(ik==0)))*om/pi*1/om1;
    g(iz)=g(iz)-u(izl)*ikl*(2-(1*(ik==0)))*om/pi;
    
    ii(end+1)=ix; jj(end+1)=ixl; vv(end+1)=-ikl*(2-(1*(ik==0)))*om/pi*1/om1;
    ii(end+1)=iy; jj(end+1)=iyl; vv(end+1)=-ikl*(2-(1*(ik==0)))*om/pi*1/om1;
    ii(end+1)=iz; jj(end+1)=izl; vv(end+1)=-ikl*(2-(1*(ik==0)))*om/pi;
    
    ii(end+1)=ix; jj(end+1)=lastid; vv(end+1)=-u(ixl)*ikl*(2-(1*(ik==0)))/pi*1/om1;
    ii(end+1)=iy; jj(end+1)=lastid; vv(end+1)=-u(iyl)*ikl*(2-(1*(ik==0)))/pi*1/om1;
    ii(end+1)=iz; jj(end+1)=lastid; vv(end+1)=-u(izl)*ikl*(2-(1*(ik==0)))/pi;
                
    end
    
    
    
end
% g(lastid)=sum(u(neq*nt+1:neq:end-1).*[0:1:nt-1]');
% aa=neq*nt+1:neq:lastid-1;
% ii=[ii,aa./aa*lastid]; jj=[jj,aa]; vv=[vv,[0:1:nt-1]];

%evp
for ik=nt-1
    
    ix=ik*neq+1; iy=ix+1; iz=iy+1;
    
    %has to express u(1)-u(0)=0
    for ikl=0:nt-1
    ixl=abs(ikl)*neq+1; iyl=abs(ikl)*neq+2; izl=abs(ikl)*neq+3; 
  
    g(ix)=g(ix)+u(ixl)-(-1)^ikl*u(ixl);
    g(iy)=g(iy)+u(iyl)-(-1)^ikl*u(iyl);
    g(iz)=g(iz)+u(izl)-(-1)^ikl*u(izl);
    
    ii(end+1)=ix; jj(end+1)=ixl; vv(end+1)=(1-(-1)^ikl);
    ii(end+1)=iy; jj(end+1)=iyl; vv(end+1)=(1-(-1)^ikl);
    ii(end+1)=iz; jj(end+1)=izl; vv(end+1)=(1-(-1)^ikl);
    
    end
end
    
    
%phase
for ikl=1:nt-1
    ixl=abs(ikl)*neq+1; %iyl=abs(ikl)*neq+2; izl=abs(ikl)*neq+3; 
  
    g(lastid)=g(lastid)+u(ixl)*ikl^2; 
    
    ii(end+1)=lastid; jj(end+1)=ixl; vv(end+1)=ikl^2; 
                
end

ii(end+1)=lastid;jj(end+1)=lastid;vv(end+1)=0;
jac=sparse(ii,jj,vv);
end