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

lastid=nt*neq*2+1;
% om=u(lastid);
for ik=0:2*nt-1
    ix=ik*neq+1; iy=ix+1; iz=iy+1;
%     
    realEq=ik<nt;
%     
%     
%     nlxz=0;
%     nlxy=0;
    kSum=mod(ik,nt);
    for ikl=kSum-nt+1:nt-1
        ikr=kSum-ikl;
        
        ixrr=abs(ikr)*neq+1; ixri=(abs(ikr)+nt)*neq+1; ixlr=abs(ikl)*neq+1; ixli=(abs(ikl)+nt)*neq+1;
        iyrr=abs(ikr)*neq+2; iyri=(abs(ikr)+nt)*neq+2; iylr=abs(ikl)*neq+2; iyli=(abs(ikl)+nt)*neq+2;
        izrr=abs(ikr)*neq+3; izri=(abs(ikr)+nt)*neq+3; izlr=abs(ikl)*neq+3; izli=(abs(ikl)+nt)*neq+3;
        
        if(realEq==1)
%             nlxz=nlxz+u(ixrr)*u(izlr)- sign(ikr)* sign(ikl)*u(ixli)*u(izri);
            g(ix)=g(ix)+u(ixrr)*u(ixlr)- sign(ikr)* sign(ikl)*u(ixli)*u(ixri);
            %jac
            ii(end+1)=ix; jj(end+1)=ixrr; vv(end+1)=u(ixlr);
            ii(end+1)=ix; jj(end+1)=ixlr; vv(end+1)=u(ixrr);
            ii(end+1)=ix; jj(end+1)=ixli; vv(end+1)=- sign(ikr)* sign(ikl)*u(ixri);
            ii(end+1)=ix; jj(end+1)=ixri; vv(end+1)=- sign(ikr)* sign(ikl)*u(ixli);
            
           
            
        else
%             nlxz=nlxz+u(ixrr)*u(izli)* sign(ikl)+ sign(ikr)*u(ixri)*u(izlr);
            g(ix)=g(ix)+u(ixrr)*u(ixli)* sign(ikl)+ sign(ikr)*u(ixri)*u(ixlr);
            
            %jac
            ii(end+1)=ix; jj(end+1)=ixrr; vv(end+1)=u(ixli)* sign(ikl);
            ii(end+1)=ix; jj(end+1)=ixli; vv(end+1)=u(ixrr)* sign(ikl);
            ii(end+1)=ix; jj(end+1)=ixri; vv(end+1)= sign(ikr)*u(ixlr);
            ii(end+1)=ix; jj(end+1)=ixlr; vv(end+1)=sign(ikr)*u(ixri);
          
        end
    end
    
    
    
    
    ikl=kSum;
    %     ixrr=abs(ikr)*neq+1; ixri=(abs(ikr)+nt)*neq+1;
    ixlr=abs(ikl)*neq+1; ixli=(abs(ikl)+nt)*neq+1;
    %         iyrr=abs(ikr)*neq+2; iyri=(abs(ikr)+nt)*neq+2;
    iylr=abs(ikl)*neq+2; iyli=(abs(ikl)+nt)*neq+2;
    %     izrr=abs(ikr)*neq+3; izri=(abs(ikr)+nt)*neq+3;
    izlr=abs(ikl)*neq+3; izli=(abs(ikl)+nt)*neq+3;
    
%     g(ix)=g(ix)+u(ix)-1*(ik==0)-0.5*(ik==1);
    g(ix)=g(ix)-3/2*(ik==0)-1*(ik==1)-0.25*(ik==2)+u(ix)*(ik==nt);
%     g(iy)=r*u(ix)-u(iy)-nlxz-(realEq*(-om*kSum*u(iyli))+~realEq*(om*kSum*u(iylr)));
%     g(iz)=nlxy-b*u(iz)-(realEq*(-om*kSum*u(izli))+~realEq*(om*kSum*u(izlr)));

   
    
    %jac
    ii(end+1)=ix; jj(end+1)=ix; vv(end+1)=(ik==nt);
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
jac=sparse(ii,jj,vv);
end