function [g,jac]=calculateRhsAndJac(neq,nt,u,r)
%%

s=10;
b=8/3;
% r=24;
ii=[]; jj=[]; vv=[];

g=u*0;

%u=[xr0,yr0,zr0,xr1,...,xi0,yi0,zi0,...]

%   dxdt=s*(y-x)
%   dydt=x*(r-z)-y
%   dzdt=xy-b*z

lastid=nt*neq+1;
om=u(lastid);
for ik=0:nt-1%2*nt-1
    ix=ik*neq+1; iy=ix+1; iz=iy+1;
    
%     nlxz=0;
%     nlxy=0;
    
    
%linear terms
    g(ix)=g(ix)+(u(iy));
    g(iy)=g(iy)-1*u(ix);
    g(iz)=g(iz)+0;
    
    ii(end+1)=ix; jj(end+1)=iy; vv(end+1)=1;
    ii(end+1)=iy; jj(end+1)=ix; vv(end+1)=-1;
    
%temporal term
    for ikl=ik+1:2:nt-1
    ixl=abs(ikl)*neq+1; iyl=abs(ikl)*neq+2; izl=abs(ikl)*neq+3; 
  
    g(ix)=g(ix)-u(ixl)*ikl*(2-(1*(ik==0)))*om/pi;
    g(iy)=g(iy)-u(iyl)*ikl*(2-(1*(ik==0)))*om/pi;
    g(iz)=g(iz)-u(izl)*ikl*(2-(1*(ik==0)))*om/pi;
    
    ii(end+1)=ix; jj(end+1)=ixl; vv(end+1)=-ikl*(2-(1*(ik==0)))*om/pi;
    ii(end+1)=iy; jj(end+1)=iyl; vv(end+1)=-ikl*(2-(1*(ik==0)))*om/pi;
    ii(end+1)=iz; jj(end+1)=izl; vv(end+1)=-ikl*(2-(1*(ik==0)))*om/pi;
    
    ii(end+1)=ix; jj(end+1)=lastid; vv(end+1)=-u(ixl)*ikl*(2-(1*(ik==0)))/pi;
    ii(end+1)=iy; jj(end+1)=lastid; vv(end+1)=-u(iyl)*ikl*(2-(1*(ik==0)))/pi;
    ii(end+1)=iz; jj(end+1)=lastid; vv(end+1)=-u(izl)*ikl*(2-(1*(ik==0)))/pi;
                
    end
    
    
%     %temp jac
%     ii(end+1)=ix; jj(end+1)=lastid; vv(end+1)=-(realEq*(-kSum*u(ixli))+~realEq*(kSum*u(ixlr)));
%     ii(end+1)=iy; jj(end+1)=lastid; vv(end+1)=-(realEq*(-kSum*u(iyli))+~realEq*(kSum*u(iylr)));
%     ii(end+1)=iz; jj(end+1)=lastid; vv(end+1)=-(realEq*(-kSum*u(izli))+~realEq*(kSum*u(izlr)));
%     
%     %jac
%     ii(end+1)=ix; jj(end+1)=ix; vv(end+1)=-s;
%     ii(end+1)=ix; jj(end+1)=iy; vv(end+1)=s;
%     ii(end+1)=ix; jj(end+1)=ixli; vv(end+1)=-(realEq*(-om*kSum));
%     ii(end+1)=ix; jj(end+1)=ixlr; vv(end+1)=-(~realEq*(om*kSum));
%     
%     ii(end+1)=iy; jj(end+1)=ix; vv(end+1)=r;
%     ii(end+1)=iy; jj(end+1)=iy; vv(end+1)=-1;
%     ii(end+1)=iy; jj(end+1)=iyli; vv(end+1)=-(realEq*(-om*kSum));
%     ii(end+1)=iy; jj(end+1)=iylr; vv(end+1)=-(~realEq*(om*kSum));
%     
%     ii(end+1)=iz; jj(end+1)=iz; vv(end+1)=-b;
%     ii(end+1)=iz; jj(end+1)=izli; vv(end+1)= -(realEq*(-om*kSum));
%     ii(end+1)=iz; jj(end+1)=izlr; vv(end+1)= -(~realEq*(om*kSum));
    
    
    
end
% g(lastid)=sum(u(neq*nt+1:neq:end-1).*[0:1:nt-1]');
% aa=neq*nt+1:neq:lastid-1;
% ii=[ii,aa./aa*lastid]; jj=[jj,aa]; vv=[vv,[0:1:nt-1]];

%phase
for ikl=1:nt-1
    ixl=abs(ikl)*neq+1; %iyl=abs(ikl)*neq+2; izl=abs(ikl)*neq+3; 
  
    g(lastid)=g(lastid)+u(ixl)*ikl^2; 
    
    ii(end+1)=lastid; jj(end+1)=ixl; vv(end+1)=ikl^2; 
                
end

ii(end+1)=lastid;jj(end+1)=lastid;vv(end+1)=0;
jac=sparse(ii,jj,vv);
end